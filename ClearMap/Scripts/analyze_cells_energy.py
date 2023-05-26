# Script with utility functions for cell count and energy calculation
import os
import numpy as np
import pandas as pd
import yaml
from yaml import Loader
from scipy.stats import ttest_ind, mannwhitneyu
import utils_PLS as upls

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

def clean_volumes_database(volumes):
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
    volumes['mean_volume'].replace(to_replace='#N/D', value=np.nan, inplace=True)
    # force values to be floats
    volumes["mean_volume"] = pd.to_numeric(volumes["mean_volume"])
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
#     # remove areas that are in macroareas 'Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei'
#     df_levels = upls.create_df_levels(volumes)
#     macroareas_to_remove = ['Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei']
#     list_areas_to_keep = df_levels[~df_levels['name_parent_l5'].isin(macroareas_to_remove)]['name_area'].values
#     df = df[df['area_name'].isin(list_areas_to_keep)]
    return df


def calculate_cells_per_area(df_mouse, area):
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


def aggregate_cells_per_area(df_mouse, vol, area):
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
        n_cells += calculate_cells_per_area(df_mouse=df_mouse, area=child)
    return n_cells


def calculate_energy_per_area(df_mouse, area, vol):
    energy = df_mouse.loc[lambda df_mouse: df_mouse['area_name'] == area, :]['source'].sum()/ \
             vol.loc[lambda vol: vol['safe_name'] == area, :]['mean_volume']
    return energy


def aggregate_energy_per_area(df_mouse, vol, area):
    # find area index
    area_id = vol.loc[lambda vol: vol['safe_name'] == area, :]['id'].values[0]
    # find area depth (IMPORTANT : not st_level!)
    area_depth = vol.loc[lambda vol: vol['safe_name'] == area, :]['depth'].values[0]
    children = find_children(area_id=area_id, l=area_depth, vol=vol)
    # loop over children
    # the first child is the area itself, and its energy is the number of neurons over the volume
    area_energy = calculate_energy_per_area(df_mouse=df_mouse, area=children[0], vol=vol).values[0]
    for child in children[1:]:
        # the other children count just by the sum of the intensities over the volume of the mother area
        area_energy += (df_mouse.loc[lambda df_mouse: df_mouse['area_name'] == child, :]['source'].sum()/ \
             vol.loc[lambda vol: vol['safe_name'] == area, :]['mean_volume']).values[0]
    return area_energy


def calculate_value_across_groups(experimental_groups, dict_results_across_mice, value='n_cells'):
    """
    Value can either be n_cells or energy
    """
    df_control = pd.DataFrame()
    df_fam = pd.DataFrame()
    df_unfam = pd.DataFrame()
    for subject in experimental_groups['Control']:
        df_control['area'] = dict_results_across_mice[subject]['area']
        df_control[subject] = dict_results_across_mice[subject][value]
    for subject in experimental_groups['Fam']:
        df_fam['area'] = dict_results_across_mice[subject]['area']
        df_fam[subject] = dict_results_across_mice[subject][value]
    for subject in experimental_groups['Unfam']:
        df_unfam['area'] = dict_results_across_mice[subject]['area']
        df_unfam[subject] = dict_results_across_mice[subject][value]
    return df_control, df_fam, df_unfam


def calculate_cells_energy_per_level(df_mouse, vol, level):
    # find all areas of a given level
    area_list = vol.loc[lambda vol: vol['st_level'] == level]['safe_name'].values
    # loop over all areas of given evel and aggregate cells per area
    n_cells_list = [aggregate_cells_per_area(df_mouse, vol, area) for area in area_list]
    # loop over all areas of given level and aggregate cells per area
    energy_list = [aggregate_energy_per_area(df_mouse=df_mouse, vol=vol, area=area) for area in area_list]
    # return df with area name and n_cells
    df_cells_energy = {'area': area_list,
                       'n_cells': n_cells_list,
                       'energy': energy_list
                       }
    df_cells_energy= pd.DataFrame(df_cells_energy)
    # remove lowercase areas
    df_cells_energy = df_cells_energy[df_cells_energy.area.str[0].str.isupper()]
    # remove areas that are in macroareas 'Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei'
    df_levels = upls.create_df_levels(vol)
    macroareas_to_remove = ['Pons', 'Medulla', 'Cerebellar cortex', 'Cerebellar nuclei']
    list_areas_to_keep = df_levels[~df_levels['name_parent_l5'].isin(macroareas_to_remove)]['name_area'].values
    df_cells_energy = df_cells_energy[df_cells_energy['area'].isin(list_areas_to_keep)]
    return df_cells_energy


def test_across_groups(df_control, df_fam, df_unfam, test='ttest'):
    """
    Test can either be ttest or mannwhitneyu
    """
    df_test = pd.DataFrame(columns=['area', 'pval_Control_vs_Fam', 
                                     'pval_Control_vs_Unfam', 'pval_Fam_vs_Unfam'])
    df_test['area'] = df_control['area']
    # loop over areas
    for area in df_control['area'].values:
        # compare control and fam
        # check if all values are not zero
        if df_control[df_control['area'] == area].values[0][1:].sum() + \
        df_fam[df_fam['area'] == area].values[0][1:].sum() != 0:
            if test == 'ttest':
                pval_control_fam = ttest_ind(df_control[df_control['area'] == area].values[0][1:],
                                             df_fam[df_fam['area'] == area].values[0][1:])
            elif test == 'mannwhitneyu':
                pval_control_fam = mannwhitneyu(df_control[df_control['area'] == area].values[0][1:],
                                                df_fam[df_fam['area'] == area].values[0][1:])
            else:
                raise ValueError('Test can either be ttest or mannwhitneyu')
            # assign pvalue to dataframe
            df_test['pval_Control_vs_Fam'][df_test.loc[df_test['area'] == area].index[0]] = pval_control_fam[1]

        # compare control and unfam
        # check if all values are not zero
        if df_control[df_control['area'] == area].values[0][1:].sum() + \
        df_unfam[df_unfam['area'] == area].values[0][1:].sum() != 0:
            if test == 'ttest':
                pval_control_unfam = ttest_ind(df_control[df_control['area'] == area].values[0][1:],
                     df_unfam[df_unfam['area'] == area].values[0][1:])
            else:
                pval_control_unfam = mannwhitneyu(df_control[df_control['area'] == area].values[0][1:],
                     df_unfam[df_unfam['area'] == area].values[0][1:])
            # assign pvalue to dataframe
            df_test['pval_Control_vs_Unfam'][df_test.loc[df_test['area'] == area].index[0]] = pval_control_unfam[1]

        # compare fam and unfam
        # check if all values are not zero
        if df_fam[df_fam['area'] == area].values[0][1:].sum() + \
        df_unfam[df_unfam['area'] == area].values[0][1:].sum() != 0:
            if test =='ttest':
                pval_fam_unfam = ttest_ind(df_fam[df_fam['area'] == area].values[0][1:],
                     df_unfam[df_unfam['area'] == area].values[0][1:])
            else:
                pval_fam_unfam = mannwhitneyu(df_fam[df_fam['area'] == area].values[0][1:],
                     df_unfam[df_unfam['area'] == area].values[0][1:])
            # assign pvalue to dataframe
            df_test['pval_Fam_vs_Unfam'][df_test.loc[df_test['area'] == area].index[0]] = pval_fam_unfam[1]
    return df_test


def cross_corr(df):
    # remove areas where no cells have been detected in any mouse
    # and remove rows with all nans
    corr_matrix = df.set_index('area').loc[
        ~(df.set_index('area')==0).all(axis=1)].dropna(axis=0).T.corr(method='pearson')
    return corr_matrix
