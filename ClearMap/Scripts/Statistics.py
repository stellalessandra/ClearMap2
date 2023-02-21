#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:56:36 2021

@author: user
"""


#%% Initialize workspace
  
  from ClearMap.Environment import *  #analysis:ignore
  import ClearMap.Analysis.Statistics.GroupStatistics as gs
  
  
  group1 = ['/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F1Control/density_counts.tif'];          
            
  group2 = ['/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F2Unfam/density_counts.tif'];
  
  
#First run the statistics over voxelized data for volumetric comparison          
  g1=gs.read_group(group1)
  g2=gs.read_group(group2)          
            
  pvals, sign=gs.t_test_voxelization(group1, group2, signed = True, remove_nan = True, p_cutoff = 0.05)
  
  pvc=gs.color_p_values(pvals, sign, positive = [1,0], negative = [0,1], p_cutoff = None, positive_trend = None, negative_trend = None, pmax = None)
  
  io.write('/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F1Control/VoxelStatisticsNeg.tif',pvc[:,:,:,0]) #Group2<Group1
  io.write('/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F2Unfam/VoxelStatisticsPos.tif',pvc[:,:,:,1]) #Group2>Group1
  
  g1avg=gs.mean(g1)
  g2avg=gs.mean(g2)
  
  io.write('/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F1Control/AVG_Control.tif',g1avg)
  io.write('/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F2Unfam/AVG_Exposed.tif',g2avg)
  
 #%%
#Run statistics over brain areas to find the ones with a significant difference  
  
  group1 = ['/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F1Control/cells_ClearMap1_700_points_transformed.npy'];
            
  group2 = ['/data01/astella/projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F2Unfam/cells_ClearMap1_700_points_transformed.npy'];
            
  points2=gs.read_group(group2,combine=False)

  points1=gs.read_group(group1,combine=False)
  
  #Set the proper annotation file
  annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(1,256)), orientation=(1,2,3),
      overwrite=False, verbose=True);
  
  counts2=gs.count_points_group_in_regions(points2, annotation_file = annotation_file, weight_group = None, invalid = 0, hierarchical = False)
  counts1=gs.count_points_group_in_regions(points1, annotation_file = annotation_file, weight_group = None, invalid = 0, hierarchical = False)
  
  labels=ano.annotation.names
  
  pvalreg, signreg=gs.t_test_region_countss(counts1, counts2, annotation_file = annotation_file, signed = True, remove_nan = True, p_cutoff = 0.05, equal_var = False)
  
  
        
