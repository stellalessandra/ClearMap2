#%% Initialize workspace
import sys
sys.path.append('/data/user/ClearMap2')
sys.path.append('/data01/astella/ClearMap2')
from ClearMap.Environment import *  #analysis:ignore
import ClearMap.Analysis.Statistics.GroupStatistics as gs
import yaml
from yaml import Loader
import os
import ClearMap.Scripts.utils as utils
import itertools

with open("ClearMap/Scripts/configfile.yaml", 'r') as stream:
  config = yaml.load(stream, Loader=Loader)

directory = config['data_directory']
subjects = [name for name in os.listdir(directory) \
              if os.path.isdir(os.path.join(directory, name))]
print(subjects)
# remove subjects
subjects.remove('BL55Control')
subjects.remove('BL76Bedding')

#%% Generate string pairs

def generate_unique_string_pairs(string_list):
    # Get all unique combinations of pairs (2 elements) from the input list
    pairs = list(itertools.combinations(string_list, 2))
    
    # Convert each pair to a tuple and remove duplicates
    unique_pairs = set(tuple(sorted(pair)) for pair in pairs)
    
    return unique_pairs
    
#%% T-test function

def t_test(group1, group2, label1, label2, directory): 
    #First run the statistics over voxelized data for volumetric comparison 
    g1 = gs.read_group(group1)
    g2 = gs.read_group(group2)
    pvals, sign=gs.t_test_voxelization(group1, group2, signed = True, 
                                        remove_nan = True, p_cutoff = 0.05)
    
    pvc=gs.color_p_values(pvals, sign, positive = [1,0], negative = [0,1], \
                          p_cutoff = None, positive_trend = None, \
                          negative_trend = None, pmax = None)
    
    g1avg=gs.mean(g1)
    g2avg=gs.mean(g2)
    
    io.write(directory+'VoxelStatisticsNeg-'+label1+'_'+label2+'.tif',pvc[:,:,:,0]) #Group2<Group1
    io.write(directory+'VoxelStatisticsPos-'+label1+'_'+label2+'.tif',pvc[:,:,:,1]) #Group2>Group1
    
    io.write(directory+'AVG_'+label1+'.tif',g1avg)
    io.write(directory+'AVG_'+label2+'.tif',g2avg)

#%% T-test Different groups
    
experimental_groups = utils.divide_in_exp_groups(list_subjects=subjects,
                                                 group_labels=['Control',
                                                               'Bedding',
                                                               'USVC57',
                                                               'USVBALB',
                                                               'USVC57_Bedding',
                                                               'USVBALB_Bedding'])
groups = list(experimental_groups.keys())

dfs = {group:[directory + subject + '/' + 'density_counts.tif' \
              for subject in subjects if group in utils._split_string(input_string=subject)] for group in groups}  

#%% T-test Different groups
    # print(dfs[pair[0]], dfs[pair[1]])
    # print(pair[0], pair[1])
for pair in list(generate_unique_string_pairs(groups)):
    print(pair[0], pair[1])
    
    t_test(group1=dfs[pair[0]], group2=dfs[pair[1]], label1=pair[0], 
               label2=pair[1], directory=directory)


#%%
#Run statistics over brain areas to find the ones with a significant difference

# control = [directory + subject + '/' + 'cells_ClearMap1_points_transformed.npy' for subject in subjects if 'Control' in subject] 

# fam = [directory + subject + '/' + 'cells_ClearMap1_points_transformed.npy' for subject in subjects if 'Fam' in subject]     

# unfam = [directory + subject + '/' + 'cells_ClearMap1_points_transformed.npy' for subject in subjects if 'Unfam' in subject]  


# points_control=gs.read_group(control, combine=False)
# points_fam=gs.read_group(fam, combine=False)
# points_unfam=gs.read_group(unfam, combine=False)

# def test_brain_areas(points1, points2):  
#     #Set the proper annotation file
#     annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
#     slicing=(slice(None),slice(None),slice(1,256)), orientation=(1,2,3),
#     overwrite=False, verbose=True)

#     counts1=gs.count_points_group_in_regions(points1, annotation_file = annotation_file, 
#                                           weight_group = None, invalid = 0, hierarchical = False)
#     counts2=gs.count_points_group_in_regions(points2, annotation_file = annotation_file, 
#                                           weight_group = None, invalid = 0, hierarchical = False)  
#     labels=ano.annotation.names

#     pvalreg, signreg=gs.t_test_region_countss(counts1, counts2, 
#                                               annotation_file = annotation_file, 
#                                               signed = True, remove_nan = True,
#                                               p_cutoff = 0.05, equal_var = False)
#     return pvalreg, signreg, labels, counts1, counts2
    
# pvalreg_control_fam, signreg_control_fam, labels, counts1cf, counts2cf = test_brain_areas(points1=points_control, points2=points_fam)  
# pvalreg_unfam_fam, signreg_unfam_fam, labels, counts1uf, counts2uf = test_brain_areas(points1=points_unfam, points2=points_fam)  
# pvalreg_control_unfam, signreg_control_unfam, labels, counts1cu, counts2cu = test_brain_areas(points1=points_control, points2=points_unfam)  



      
