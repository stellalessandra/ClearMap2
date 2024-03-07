# Script with utility functions for cell count and energy calculation
import os
import numpy as np
import itertools
import pandas as pd
import yaml
from yaml import Loader
from scipy.stats import ttest_ind, mannwhitneyu
import utils_PLS as upls
import scipy

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
    query = pd.read_csv("query.csv")
    table_s4 = pd.read_excel('table_S4_wang2020.xlsx', header=1)
    temp = pd.concat([table_s4['structure ID'], table_s4['CCF Volume']], axis=1)
    temp.columns = ['id', 'CCF_volume']
    volumes = pd.merge(query, temp, how='left', left_on='id', right_on='id')
    # now we split the path of a single area into columns with all the upper layers
    structure_path = volumes['structure_id_path']
    #split volume into depths
    structure_path = volumes.structure_id_path.str.split("/", expand = True)
    #remove useless columns
    structure_path=structure_path.drop(structure_path.columns[[0,12]], axis=1)
    #append data frame to volumes
    volumes = pd.concat([volumes, structure_path], axis=1)
    # rename columns to have the correct depth level
    rename_cols = {i+1:i for i in range(0,11)}
    volumes.rename(columns=rename_cols, inplace=True)
    for i in range(0,11):
        volumes[i] = pd.to_numeric(volumes[i])
    # remove a few useless columns
    columns_to_drop = ['ontology_id', 'hemisphere_id', 'weight', 'graph_id', 'graph_order', 'color_hex_triplet', \
    'neuro_name_structure_id', 'neuro_name_structure_id_path', 'failed', 'sphinx_id', 'structure_name_facet', \
     'failed_facet']
    volumes = volumes.drop(columns_to_drop, axis=1)
    # replace N/D values from excel by np.nan
    volumes['CCF_volume'].replace(to_replace='#N/D', value=np.nan, inplace=True)
    # force values to be floats
    volumes["CCF_volume"] = pd.to_numeric(volumes["CCF_volume"])
    volumes.columns
    return volumes

def reformat_df_mouse(df):
    df.columns = ['x_resampled', 'y_resampled', 'z_resampled',
                 'size', 'source', 'x', 'y', 'z', 'id', 
                  'area_name', 'layer', 'layer1', 'layer2']
    # remove slashes in the area section
    df.replace(to_replace= r'area/Layer', value= 'area Layer', regex=True, inplace=True)
    # join area names with layers
    cols = ['area_name', 'layer', 'layer1', 'layer2']
#     df['area_name'] = df[cols].apply(lambda row: ' '.join(row.values.astype(str)), axis=1)
    df['area_name'] = df['area_name'].astype(str) + df['layer'].fillna('') + ' ' \
                + df['layer1'].fillna('') + ' ' + df['layer2'].fillna('')
    #remove nans
    df.replace(to_replace= r'nan', value= '', regex=True, inplace=True)
    #remove double spaces
    df.replace(to_replace= r'  ', value= ' ', regex=True, inplace=True)
    #remove trailing spaces at the end of the string
    df['area_name'] = df['area_name'].str.rstrip()
    df = df.drop(['layer', 'layer1', 'layer2'], axis=1)
    return df


def calculate_cells_per_area(df_mouse, area, source, size):
    # apply second filtering on source or size
    if source != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['source'] > source, :]
    if size != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['size'] > size, :]
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
    # find area index
    area_id = vol.loc[lambda vol: vol['safe_name'] == area, :]['id'].values[0]
    # find area depth (IMPORTANT : not st_level!)
    area_depth = vol.loc[lambda vol: vol['safe_name'] == area, :]['depth'].values[0]
    # now I have to find all areas that have that area up as parent up to the area level
    children = find_children(area_id=area_id, l=area_depth, vol=vol)
    # loop over children
    n_cells = 0
    for child in children:
        # count number of cells per area
        n_cells += calculate_cells_per_area(df_mouse=df_mouse, area=child, source=source, size=size)
    return n_cells


def calculate_energy_per_area(df_mouse, area, vol, source, size):
    # apply second filtering on source or size
    if source != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['source'] > source, :]
    if size != 0:
        df_mouse = df_mouse.loc[lambda df_mouse: df_mouse['size'] > size, :]
    energy = df_mouse.loc[lambda df_mouse: df_mouse['area_name'] == area, :]['source'].sum()/ \
             vol.loc[lambda vol: vol['safe_name'] == area, :]['CCF_volume']
    return energy


def aggregate_energy_per_area(df_mouse, vol, area, source, size):
    # find area index
    area_id = vol.loc[lambda vol: vol['safe_name'] == area, :]['id'].values[0]
    # find area depth (IMPORTANT : not st_level!)
    area_depth = vol.loc[lambda vol: vol['safe_name'] == area, :]['depth'].values[0]
    children = find_children(area_id=area_id, l=area_depth, vol=vol)
    # loop over children
    # the first child is the area itself, and its energy is the number of neurons over the volume
    area_energy = calculate_energy_per_area(df_mouse=df_mouse, area=children[0], vol=vol, source=source, size=size).values[0]
    for child in children[1:]:
        # the other children count just by the sum of the intensities over the volume of the mother area
        area_energy += (df_mouse.loc[lambda df_mouse: df_mouse['area_name'] == child, :]['source'].sum()/ \
             vol.loc[lambda vol: vol['safe_name'] == area, :]['CCF_volume']).values[0]
    return area_energy


def calculate_value_across_groups(experimental_groups, dict_results_across_mice, value='n_cells'):
    """
    Value can either be n_cells or energy
    """
    dfs = [pd.DataFrame() for i in range(len(experimental_groups.keys()))]
    
    for i, group in enumerate(experimental_groups.keys()):
        for subject in experimental_groups[group]:
            dfs[i]['area'] = dict_results_across_mice[subject]['area']
            dfs[i][subject] = dict_results_across_mice[subject][value]
    return dfs


def calculate_cells_energy_per_level(df_mouse, vol, level, source=0, size=0):
    # find all areas of a given level
    area_list = vol.loc[lambda vol: vol['st_level'] == level]['safe_name'].values
    # loop over all areas of given evel and aggregate cells per area
    n_cells_list = [aggregate_cells_per_area(df_mouse=df_mouse, vol=vol, area=area, source=source, size=size) for area in area_list]
    # loop over all areas of given level and aggregate cells per area
    energy_list = [aggregate_energy_per_area(df_mouse=df_mouse, vol=vol, area=area, source=source, size=size) for area in area_list]
    # return df with area name and n_cells
    df = {'area': area_list,
                       'n_cells': n_cells_list,
                       'energy': energy_list
                       }
    df= pd.DataFrame(df)
    # remove lowercase areas
    df = df[df.area.str[0].str.isupper()]
    # remove areas that are in macroareas 'Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei'
    df_levels = upls.create_df_levels(vol)
    macroareas_to_remove = ['Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei']
    list_areas_to_keep = df_levels[~df_levels['name_parent_l5'].isin(macroareas_to_remove)]['name_area'].values
    df = df[df['area'].isin(list_areas_to_keep)]
    
    #loop over all areas of given level and calculate density
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
    return df


def test_across_groups(dfs, groups=['Control', 'Fam', 'Unfam'], test='ttest'):
    """
    Test can either be ttest or mannwhitneyu
    """
    #calculate groups combinations and column names for dataframe
    combinations = list(itertools.combinations(groups, 2))
    pval_keys = ['pval_'+pair[0]+'_vs_'+pair[1] for pair in combinations]
    df_test = pd.DataFrame(columns=['area'] + pval_keys)
    #fix area list
    df_test['area'] = dfs[0]['area']
    # loop over areas
    for area in dfs[0]['area'].values:
        for i, (df1, df2) in enumerate(list(itertools.combinations(dfs, 2))):
            # compare control and fam
            # check if all values are not zero
            if df1[df1['area'] == area].values[0][1:].sum() + \
            df2[df2['area'] == area].values[0][1:].sum() != 0:
                if test == 'ttest':
                    pval = ttest_ind(df1[df1['area'] == area].values[0][1:],
                                     df2[df2['area'] == area].values[0][1:])
                elif test == 'mannwhitneyu':
                    pval = mannwhitneyu(df1[df1['area'] == area].values[0][1:],
                                        df2[df2['area'] == area].values[0][1:])
                else:
                    raise ValueError('Test can either be ttest or mannwhitneyu')
                # assign pvalue to dataframe
                df_test[pval_keys[i]][df_test.loc[df_test['area'] == area].index[0]] = pval[1]
    return df_test


def cross_corr(df):
    # remove areas where no cells have been detected in any mouse
    # and remove rows with all nans
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
    combinations = list(itertools.combinations(groups, 2))
    tests = ['pval_'+pair[0]+'_vs_'+pair[1] for pair in combinations]
    d = []
    for key in tests:
        d.append(df.sort_values(by=key)[[
        'area', key]]['area'][df.sort_values(by=key)[[
        'area', key]][key] < threshold].to_list())
    # return flattened list
    return np.unique(list(itertools.chain.from_iterable(d)))


def select_significant_areas(dictionary, experimental_groups, batch,  
                             test='mannwhitneyu', threshold_test=0.05,
                             threshold_pls=2.56,
                            value_test='n_cells', value_pls='relative_density'):
    
    volumes = clean_volumes_database()
    dfs = \
    calculate_value_across_groups(experimental_groups=experimental_groups, 
                                  dict_results_across_mice=dictionary, 
                                  value=value_test)
    groups = list(experimental_groups.keys())
    df_test = test_across_groups(dfs,
                                 test=test,
                                groups=groups)

    combinations = list(itertools.combinations(groups, 2))
    tests = ['pval_'+pair[0]+'_vs_'+pair[1] for pair in combinations]
    
    d = []
    for key in tests:
        d.append(df_test.sort_values(by=key)[[
        'area', key]]['area'][df_test.sort_values(by=key)[[
        'area', key]][key] < threshold_test].to_list())
    # return flattened list
    sig_areas_test = np.unique(list(itertools.chain.from_iterable(d)))
    
    #now find significant areas for PLS
    overlap = {'ncells':[], 'energy':[], 'density':[], 'relative_density':[]}
    for variable in overlap.keys():
        overlap[variable] = set(upls.identify_pls_sig_areas(saliences=pd.read_csv(
            './results_pls/'+ batch +'_'+ variable +'_saliences.csv'), 
                                               threshold=threshold_pls, 
                                               volumes=clean_volumes_database()))
    significant_areas = overlap[value_pls].union(list(sig_areas_test))
    return significant_areas
    
    
    
