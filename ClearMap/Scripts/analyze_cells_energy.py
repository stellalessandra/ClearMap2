import os
import numpy as np
import itertools
import pandas as pd
import yaml
from yaml import Loader
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
import utils_PLS as upls
import scipy
import scikit_posthocs as sp


def list_subjects(root_directory):
    try:
        with open("ClearMap/Scripts/configfile.yaml", 'r') as stream:
            config = yaml.load(stream, Loader=Loader)
    except:
        with open("configfile.yaml", 'r') as stream:
            config = yaml.load(stream, Loader=Loader)
    
    user = config['user']
    experiment = config['experiment']
    experimental_group = config['experimental_group']
    
    data_directory = root_directory + experiment + '/' \
                + experimental_group + '/'
                
    subjects = [name for name in os.listdir(data_directory) \
                if os.path.isdir(os.path.join(data_directory, name))]
    
    return subjects

def clean_volumes_database():
    """
    This function cleans the volumes database by reading data from a CSV and an Excel file, merging them, 
    splitting the structure path into depths, renaming columns, and removing unnecessary columns. 
    It also handles missing values and ensures that the 'CCF_volume' column is of float type.
    
    Returns:
    pd.DataFrame: A cleaned pandas DataFrame representing the volumes database.
    """
    
    # Read the query data from a CSV file
    query = pd.read_csv("query.csv")
    
    # Read the table data from an Excel file
    table_s4 = pd.read_excel('table_S4_wang2020.xlsx', header=1)
    
    # Concatenate the 'structure ID' and 'CCF Volume' columns from the table
    temp = pd.concat([table_s4['structure ID'], table_s4['CCF Volume']], axis=1)
    temp.columns = ['id', 'CCF_volume']
    
    # Merge the query and temp DataFrames based on 'id'
    volumes = pd.merge(query, temp, how='left', left_on='id', right_on='id')
    
    # Split the 'structure_id_path' into depths
    structure_path = volumes.structure_id_path.str.split("/", expand = True)
    
    # Remove the first and last columns as they are not needed
    structure_path=structure_path.drop(structure_path.columns[[0,12]], axis=1)
    
    # Append the structure_path DataFrame to the volumes DataFrame
    volumes = pd.concat([volumes, structure_path], axis=1)
    
    # Rename the columns to have the correct depth level
    rename_cols = {i+1:i for i in range(0,11)}
    volumes.rename(columns=rename_cols, inplace=True)
    
    # Convert the depth columns to numeric
    for i in range(0,11):
        volumes[i] = pd.to_numeric(volumes[i])
    
    # Remove unnecessary columns
    columns_to_drop = ['ontology_id', 'hemisphere_id', 'weight', 'graph_id', 'graph_order', 'color_hex_triplet', \
    'neuro_name_structure_id', 'neuro_name_structure_id_path', 'failed', 'sphinx_id', 'structure_name_facet', \
     'failed_facet']
    volumes = volumes.drop(columns_to_drop, axis=1)
    
    # Replace '#N/D' values in the 'CCF_volume' column with np.nan
    volumes['CCF_volume'].replace(to_replace='#N/D', value=np.nan, inplace=True)
    
    # Convert the 'CCF_volume' column to float
    volumes["CCF_volume"] = pd.to_numeric(volumes["CCF_volume"])
    
    return volumes


def reformat_df_mouse(df):
    """
    This function reformats a given DataFrame by renaming columns, 
    cleaning and merging certain fields, and removing unnecessary columns.

    Parameters:
    df (pd.DataFrame): Input DataFrame to be reformatted.

    Returns:
    df (pd.DataFrame): Reformatted DataFrame.
    """

    # Rename the columns of the DataFrame
    df.columns = ['x_resampled', 'y_resampled', 'z_resampled',
                 'size', 'source', 'x', 'y', 'z', 'id', 
                  'area_name', 'layer', 'layer1', 'layer2']

    # Remove slashes in the 'area' section
    df.replace(to_replace= r'area/Layer', value= 'area Layer', regex=True, inplace=True)

    # Join 'area_name' with 'layer', 'layer1', and 'layer2'
    # The 'fillna('')' function is used to replace NaN values with an empty string before joining
    df['area_name'] = df['area_name'].astype(str) + df['layer'].fillna('') + ' ' \
                + df['layer1'].fillna('') + ' ' + df['layer2'].fillna('')

    # Remove 'nan' strings resulting from the previous operation
    df.replace(to_replace= r'nan', value= '', regex=True, inplace=True)

    # Remove double spaces that might have been introduced in the 'area_name' field
    df.replace(to_replace= r'  ', value= ' ', regex=True, inplace=True)

    # Remove trailing spaces at the end of the 'area_name' string
    df['area_name'] = df['area_name'].str.rstrip()

    # Drop the now unnecessary 'layer', 'layer1', and 'layer2' columns
    df = df.drop(['layer', 'layer1', 'layer2'], axis=1)

    return df


def calculate_cells_per_area(df_mouse, area, source, size):
    """
    This function calculates the number of cells in a specific area after applying filters on 'source' and 'size'.
    
    Parameters:
    df_mouse (DataFrame): The input DataFrame containing information about the mouse cells.
    area (str): The specific area to calculate the number of cells for.
    source (int): The minimum source value to filter the DataFrame. If it's 0, no filtering is applied.
    size (int): The minimum size value to filter the DataFrame. If it's 0, no filtering is applied.

    Returns:
    n_cells (int): The number of cells in the specified area after applying the filters.
    """
    
    # Apply second filtering on source if source is not 0
    if source != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['source'] > source, :]
    
    # Apply second filtering on size if size is not 0
    if size != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['size'] > size, :]
    
    # Calculate the number of cells in the specified area
    n_cells = len(df_mouse.loc[lambda df_mouse: df_mouse['area_name'] == area, :])
    
    return n_cells


def find_children(area_id, l, vol):
    '''
    Function that finds, given one area id, all areas having that area 
    in their hierarchy up to level l
    input: area_id, level, volume table
    '''
    children = []
    for i in range(l, 11):
        found_subareas = vol.loc[lambda vol: vol[i] == area_id, :]['safe_name'].values
        if len(found_subareas):
            children.append(found_subareas)
    children = np.array(children).flatten()
    return children



def aggregate_cells_per_area(df_mouse, vol, area, source, size):
    """
    This function aggregates the number of cells per area.

    Parameters:
    df_mouse (DataFrame): The DataFrame containing mouse data.
    vol (DataFrame): The DataFrame containing volume data.
    area (str): The name of the area.
    source (int): The source value for filtering.
    size (int): The size value for filtering.

    Returns:
    int: The total number of cells in the specified area.
    """

    # Find the index of the specified area
    area_id = vol.loc[lambda vol: vol['safe_name'] == area, :]['id'].values[0]

    # Find the depth of the specified area (IMPORTANT : not st_level!)
    area_depth = vol.loc[lambda vol: vol['safe_name'] == area, :]['depth'].values[0]

    # Find all areas that have the specified area as a parent up to the area level
    children = find_children(area_id=area_id, l=area_depth, vol=vol)

    # Initialize the cell counter
    n_cells = 0

    # Loop over each child area
    for child in children:
        # Count the number of cells in the child area and add it to the total
        n_cells += calculate_cells_per_area(df_mouse=df_mouse, area=child, source=source, size=size)

    # Return the total number of cells
    return n_cells


def calculate_energy_per_area(df_mouse, area, vol, source, size):
    """
    This function calculates the energy per area.

    Parameters:
    df_mouse (DataFrame): The DataFrame containing mouse data.
    area (str): The name of the area.
    vol (DataFrame): The DataFrame containing volume data.
    source (int): The source value for filtering.
    size (int): The size value for filtering.

    Returns:
    float: The energy per area.
    """

    # Apply a filter on the source if it is not zero
    if source != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['source'] > source, :]

    # Apply a filter on the size if it is not zero
    if size != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['size'] > size, :]

    # Calculate the energy by summing the source values for the specified area and dividing by the volume of the area
    energy = df_mouse.loc[lambda df_mouse: df_mouse['area_name'] == area, :]['source'].sum() / \
             vol.loc[lambda vol: vol['safe_name'] == area, :]['CCF_volume']

    # Return the energy per area
    return energy


def aggregate_energy_per_area(df_mouse, vol, area, source, size):
    """
    This function aggregates the energy per area.

    Parameters:
    df_mouse (DataFrame): The DataFrame containing mouse data.
    vol (DataFrame): The DataFrame containing volume data.
    area (str): The name of the area.
    source (int): The source value for filtering.
    size (int): The size value for filtering.

    Returns:
    float: The total energy in the specified area.
    """

    # Find the index of the specified area
    area_id = vol.loc[lambda vol: vol['safe_name'] == area, :]['id'].values[0]

    # Find the depth of the specified area (IMPORTANT : not st_level!)
    area_depth = vol.loc[lambda vol: vol['safe_name'] == area, :]['depth'].values[0]

    # Find all areas that have the specified area as a parent up to the area level
    children = find_children(area_id=area_id, l=area_depth, vol=vol)

    # Calculate the energy for the first child area (the area itself)
    area_energy = calculate_energy_per_area(df_mouse=df_mouse, area=children[0], vol=vol, source=source, size=size).values[0]

    # Loop over each child area (excluding the first one)
    for child in children[1:]:
        # Add the energy of the child area to the total energy
        # The energy of a child area is calculated as the sum of the intensities over the volume of the parent area
        area_energy += (df_mouse.loc[lambda df_mouse: df_mouse['area_name'] == child, :]['source'].sum() / \
             vol.loc[lambda vol: vol['safe_name'] == area, :]['CCF_volume']).values[0]

    # Return the total energy
    return area_energy


def calculate_value_across_groups(experimental_groups, dict_results_across_mice, value='n_cells'):
    """
    This function calculates a specified value across different experimental groups.

    Parameters:
    experimental_groups (dict): A dictionary where keys are group names and values are lists of subjects in each group.
    dict_results_across_mice (dict): A dictionary where keys are subject names and values are 
    DataFrames containing results for each subject.
    value (str): The value to calculate. Can be either 'n_cells', 'energy', 'density' or 'relative density'. 
    Default is 'n_cells'.

    Returns:
    list: A list of DataFrames, one for each experimental group. Each DataFrame contains the calculated v
    alue for each subject in the group.
    """

    # Initialize a list of empty DataFrames, one for each experimental group
    dfs = [pd.DataFrame() for i in range(len(experimental_groups.keys()))]
    
    # Loop over each experimental group
    for i, group in enumerate(experimental_groups.keys()):
        # Loop over each subject in the current group
        for subject in experimental_groups[group]:
            # Add the area and calculated value for the current subject to the DataFrame for the current group
            dfs[i]['area'] = dict_results_across_mice[subject]['area']
            dfs[i][subject] = dict_results_across_mice[subject][value]

    # Return the list of DataFrames
    return dfs



def calculate_cells_energy_per_level(df_mouse, vol, level, source=0, size=0):
    """
    This function calculates the number of cells, energy, density 
    and relative density per level.

    Parameters:
    df_mouse (DataFrame): The DataFrame containing mouse data.
    vol (DataFrame): The DataFrame containing volume data.
    level (int): The level to calculate values for.
    source (int): The source value for filtering. Default is 0.
    size (int): The size value for filtering. Default is 0.

    Returns:
    DataFrame: A DataFrame containing the number of cells, energy, density, and 
    relative density for each area at the specified level.
    """

    # Find all areas at the specified level
    area_list = vol.loc[lambda vol: vol['st_level'] == level]['safe_name'].values

    # Loop over all areas at the specified level and aggregate cells per area
    n_cells_list = [aggregate_cells_per_area(df_mouse=df_mouse, vol=vol, area=area, 
                                             source=source, size=size) for area in area_list]

    # Loop over all areas at the specified level and aggregate energy per area
    energy_list = [aggregate_energy_per_area(df_mouse=df_mouse, vol=vol, area=area, 
                                             source=source, size=size) for area in area_list]

    # Create a DataFrame with area name, number of cells, and energy
    df = {'area': area_list,
          'n_cells': n_cells_list,
          'energy': energy_list
          }
    df = pd.DataFrame(df)

    # Remove lowercase areas
    df = df[df.area.str[0].str.isupper()]

    # Remove areas that are in macroareas 'Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei'
    df_levels = upls.create_df_levels(vol)
    macroareas_to_remove = ['Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei']
    list_areas_to_keep = df_levels[~df_levels['name_parent_l5'].isin(macroareas_to_remove)]['name_area'].values
    df = df[df['area'].isin(list_areas_to_keep)]
    
    # Loop over all areas at the specified level and calculate density and relative density
    total_volume = vol[vol['safe_name']=='root']['CCF_volume'].values[0]
    density_list = []
    rel_density_list = []
    for idx, area in enumerate(list_areas_to_keep):
        volume = vol.loc[lambda vol: vol['safe_name'] == area, :]['CCF_volume'].values[0]
        n_cells = df.loc[df['area'] == area]['n_cells'].values[0]
        density_list.append(n_cells/volume)
        rel_density_list.append((n_cells/volume)/(df['n_cells'].sum()/total_volume))
    df['density'] = density_list
    df['relative_density'] = rel_density_list

    # Return the DataFrame
    return df



def test_across_groups(dfs, groups=['Control', 'Fam', 'Unfam'], test='ttest'):
    """
    This function performs statistical tests across different groups.

    Parameters:
    dfs (list): A list of pandas dataframes, each representing a group.
    groups (list): A list of group names. Default is ['Control', 'Fam', 'Unfam'].
    test (str): The type of statistical test to perform. Can be 'ttest' or 'mannwhitneyu'. Default is 'ttest'.

    Returns:
    df_test (DataFrame): A pandas dataframe with p-values of the tests.
    """

    # Calculate combinations of groups and column names for the result dataframe
    combinations = list(itertools.combinations(groups, 2))
    pval_keys = ['pval_'+pair[0]+'_vs_'+pair[1] for pair in combinations]

    # Initialize the result dataframe with columns for each area and p-value
    df_test = pd.DataFrame(columns=['area'] + pval_keys)

    # Set the area list from the first dataframe in dfs
    df_test['area'] = dfs[0]['area']

    # Loop over areas
    for area in dfs[0]['area'].values:
        # Loop over combinations of dataframes (groups)
        for i, (df1, df2) in enumerate(list(itertools.combinations(dfs, 2))):
            # Check if all values are not zero for the current area in both dataframes
            if df1[df1['area'] == area].values[0][1:].sum() + \
            df2[df2['area'] == area].values[0][1:].sum() != 0:
                # Perform the specified test
                if test == 'ttest':
                    pval = ttest_ind(df1[df1['area'] == area].values[0][1:],
                                     df2[df2['area'] == area].values[0][1:])
                elif test == 'mannwhitneyu':
                    pval = mannwhitneyu(df1[df1['area'] == area].values[0][1:],
                                        df2[df2['area'] == area].values[0][1:])
                else:
                    raise ValueError('Test can either be ttest or mannwhitneyu')
                # Assign the p-value to the result dataframe
                df_test[pval_keys[i]][df_test.loc[df_test['area'] == area].index[0]] = pval[1]
    return df_test


def cross_corr(df):
    """
    This function calculates the Pearson correlation matrix for a given dataframe.

    Parameters:
    df (DataFrame): A pandas dataframe with 'area' as one of the columns.

    Returns:
    corr_matrix (DataFrame): A pandas dataframe representing the correlation matrix.
    """

    # Remove areas where no cells have been detected in any mouse
    # and remove rows with all NaNs
    corr_matrix = df.set_index('area').loc[
        ~(df.set_index('area')==0).all(axis=1)].dropna(axis=0).T.corr(method='pearson')
    return corr_matrix


def compute_corr_pvalue(dfs_groups, group_labels, volumes, behavioral_data, significant_areas=None,
                       correlation='Pearson'):
    """
    df_groups: list of pandas dataframe
    group_labels: list of strings
    volumes: volumes dataframe of Allen Brain Atlas
    behavioral_data: behavioral variable registered
    
    Returns
    -------
    corr_df, pvalue_df: pandas dataframes of pearson correlation and corresponding pvalue
    between the given dataframe rows and the behavioral variable
    """
    areas = [volumes[volumes['safe_name']==area]['acronym'].values[0] \
             for area in dfs_groups[0].dropna().reset_index(drop=True)['area']]

    corr_df = pd.DataFrame(columns=areas, index=group_labels)
    pvalue_df = pd.DataFrame(columns=areas, index=group_labels)
    
    for i, (df_group, group_label) in enumerate(zip(dfs_groups, group_labels)):
        corr = []
        pvalue = []
        df = df_group.dropna().reset_index(drop=True)

        for area in areas:
            area_name = volumes[volumes['acronym']==area]['safe_name'].values[0]
            if correlation=='Pearson':
                corr.append(scipy.stats.pearsonr(df.set_index('area').loc[area_name].values, 
                                 behavioral_data[i].values)[0])
                pvalue.append(scipy.stats.pearsonr(df.set_index('area').loc[area_name].values, 
                                 behavioral_data[i].values)[1])
            elif correlation=='Spearmann':
                corr.append(scipy.stats.spearmanr(df.set_index('area').loc[area_name].values, 
                                 behavioral_data[i].values)[0])
                pvalue.append(scipy.stats.spearmanr(df.set_index('area').loc[area_name].values, 
                                 behavioral_data[i].values)[1])
            elif correlation=='Kendall':
                corr.append(scipy.stats.kendalltau(df.set_index('area').loc[area_name].values, 
                                 behavioral_data[i].values)[0])
                pvalue.append(scipy.stats.kendalltau(df.set_index('area').loc[area_name].values, 
                                 behavioral_data[i].values)[1])
            else:
                raise ValueError('Wrong correlation parameter')
        corr_df.loc[group_label] = corr
        pvalue_df.loc[group_label] = pvalue
    if significant_areas is not None:
        final_areas = set(corr_df.columns).intersection(set(significant_areas))
        corr_df = corr_df[final_areas]
        pvalue_df = pvalue_df[final_areas]
    corr_df = corr_df[corr_df.columns].astype(float)
    pvalue_df = pvalue_df[pvalue_df.columns].astype(float)
    return corr_df, pvalue_df


def select_areas(df, threshold=0.05, groups=['Control', 'Fam', 'Unfam']):
    """
    This function selects areas based on a given threshold.

    Parameters:
    df (DataFrame): A pandas dataframe with 'area' as one of the columns and p-values.
    threshold (float): The threshold for p-values. Default is 0.05.
    groups (list): A list of group names. Default is ['Control', 'Fam', 'Unfam'].

    Returns:
    (numpy array): A flattened list of unique areas that meet the threshold.
    """

    # Calculate combinations of groups and create a list of test names
    combinations = list(itertools.combinations(groups, 2))
    tests = ['pval_'+pair[0]+'_vs_'+pair[1] for pair in combinations]

    # Initialize a list to store areas
    d = []

    # Loop over each test
    for key in tests:
        # Append areas where p-value is less than the threshold to the list
        d.append(df.sort_values(by=key)[['area', key]]['area'][df.sort_values(by=key)[['area', key]][key] < threshold].to_list())

    # Return a flattened list of unique areas
    return np.unique(list(itertools.chain.from_iterable(d)))



def select_significant_areas(dictionary, experimental_groups, batch,  
                             test='mannwhitneyu', threshold_test=0.05,
                             threshold_pls=2.56,
                            value_test='n_cells', value_pls='relative_density'):
    """
    This function selects significant areas based on a given threshold and test.

    Parameters:
    dictionary (dict): A dictionary of results across mice.
    experimental_groups (dict): A dictionary of experimental groups.
    batch (str): The batch name.
    test (str): The type of statistical test to perform. Can be 'ttest' or 'mannwhitneyu'. Default is 'mannwhitneyu'.
    threshold_test (float): The threshold for the test's p-values. Default is 0.05.
    threshold_pls (float): The threshold for the PLS. Default is 2.56.
    value_test (str): The value for the test. Default is 'n_cells'.
    value_pls (str): The value for the PLS. Default is 'relative_density'.

    Returns:
    significant_areas (set): A set of significant areas.
    """

    # Clean volumes database
    volumes = clean_volumes_database()

    # Calculate value across groups
    dfs = calculate_value_across_groups(experimental_groups=experimental_groups, 
                                        dict_results_across_mice=dictionary, 
                                        value=value_test)

    # List of group names
    groups = list(experimental_groups.keys())

    # Perform tests across groups
    df_test = test_across_groups(dfs, test=test, groups=groups)

    # Calculate combinations of groups and create a list of test names
    combinations = list(itertools.combinations(groups, 2))
    tests = ['pval_'+pair[0]+'_vs_'+pair[1] for pair in combinations]

    # Initialize a list to store areas
    d = []

    # Loop over each test
    for key in tests:
        # Append areas where p-value is less than the threshold to the list
        d.append(df_test.sort_values(by=key)[['area', key]]['area'][df_test.sort_values(by=key)[['area', key]][key] < threshold_test].to_list())

    # Get unique areas from the list
    sig_areas_test = np.unique(list(itertools.chain.from_iterable(d)))

    # Initialize a dictionary to store overlap
    overlap = {'ncells':[], 'energy':[], 'density':[], 'relative_density':[]}

    # Find significant areas for PLS
    for variable in overlap.keys():
        overlap[variable] = set(upls.identify_pls_sig_areas(saliences=pd.read_csv(
            './results_pls/'+ batch +'_'+ variable +'_saliences.csv'), 
                                               threshold=threshold_pls, 
                                               volumes=clean_volumes_database()))

    # Get union of significant areas from the test and PLS
    significant_areas = overlap[value_pls].union(list(sig_areas_test))

    return significant_areas


def create_dfs_across_groups(dictionary_results, experimental_groups, value, area=None):
    """
    This function creates a dictionary of pandas DataFrames for each experimental group.

    Parameters:
    dictionary_results (dict): A dictionary containing the results of the experiments.
    experimental_groups (dict): A dictionary where the keys are the names of the 
    experimental groups and the values are lists of subjects in each group.
    value (str): The value to be calculated across groups.
    area (str, optional): The specific area to be considered. If None, the function will 
    sum the values across all areas.

    Returns:
    dict: A dictionary where the keys are the names of the experimental groups and the 
    values are pandas DataFrames. Each DataFrame has columns 'subject', 'value', and 'group'. 
    The 'subject' column contains the subjects in the group, the 'group' column contains the 
    name of the group, and the 'value' column contains the calculated value for each subject in the group.

    """
    # Calculate the specified value across groups for each experimental group
    dfs_value = {key: calculate_value_across_groups(experimental_groups=experimental_groups, 
                              dict_results_across_mice=dictionary_results, 
                              value=value)[i] for i, key in enumerate(experimental_groups.keys())}
    
    # Initialize a DataFrame for each experimental group
    dfs = {key: pd.DataFrame(columns=['subject', value, 'group']) for key in experimental_groups.keys()}
    
    for group, df in dfs.items():
        # Assign the subjects and group name to the DataFrame
        df['subject'] = [s for s in experimental_groups[group]]
        df['group'] = [group for s in experimental_groups[group]]
        
        if area==None:
            # If no area is specified, sum the values across all areas
            df[value] = [dfs_value[group].set_index('area').sum()[s] for s in experimental_groups[group]]
        else:
            # If an area is specified, calculate the value for that specific area
            df[value] = [dfs_value[group].set_index('area').loc[area][s] for s in experimental_groups[group]]

    # Return the dictionary of DataFrames
    return dfs



def kruskal_per_area(dictionary, value, experimental_groups):
    """
    Performs a Kruskal-Wallis H-test on a given dictionary of dataframes for each area.

    Parameters:
    dictionary (dict): A dictionary where keys are subjects and values are pandas dataframes.
    value (str): The column name in the dataframe to perform the test on.
    experimental_groups (list): A list of experimental groups.

    Returns:
    dict: A dictionary where keys are areas and values are tuples. 
    Each tuple contains the result of the Kruskal-Wallis H-test and the result of the post-hoc Dunn's test.
    
    """
    
    # Get the list of subjects from the dictionary keys
    subjects = list(dictionary.keys())

    # Create a dataframe of levels from the cleaned volumes database
    df_levels = upls.create_df_levels(clean_volumes_database())

    # Filter out certain areas and convert the remaining areas to a numpy array
    areas = df_levels[~df_levels['parent_l5'].isin(['P','MY','CBX', 'CBN'])]['name_area'].to_numpy()

    # Initialize an empty dictionary to store significant areas
    significant_areas = {}

    # Loop over each area
    for area in areas:
        # Create dataframes for the results dictionary, experimental groups, value, and area
        dfs = create_dfs_across_groups(dictionary_results=dictionary, 
                experimental_groups=experimental_groups, 
                value=value, area=area)

        # Calculate the total sum of the value column for each key in the dataframes
        tot_sum = 0
        for key in dfs.keys():
            tot_sum += dfs[key][value].to_numpy().sum()

        # If the total sum is not zero, perform a Kruskal-Wallis H-test
        if tot_sum != 0:
            # If the p-value of the test is less than 0.05, add the area to the significant_areas dictionary
            if kruskal(*[dfs[key][value].to_numpy() for key in dfs.keys()])[1] < 0.05:
                significant_areas[area] = \
                    (kruskal(*[dfs[key][value].to_numpy() for key in dfs.keys()]), 
                     sp.posthoc_dunn(pd.concat(dfs, ignore_index=True), val_col=value, 
                group_col='group', p_adjust = 'fdr_bh'))

    return significant_areas
    
