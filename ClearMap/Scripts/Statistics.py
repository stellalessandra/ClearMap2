#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 10:56:36 2021

@author: user
"""


#%% Initialize workspace
  
  from ClearMap.Environment import *  #analysis:ignore
  import ClearMap.Analysis.Statistics.GroupStatistics as gs
  
  
  control = ['/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F1Control/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F10Control/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F11Control/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F16Control/density_counts.tif']          
            
  fam = ['/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F3Fam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F5Fam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F7Fam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F9Fam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F13Fam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F15Fam/density_counts.tif']
  
  unfam = ['/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F2Unfam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F4Unfam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F6Unfam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F8Unfam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F12Unfam/density_counts.tif',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F14Unfam/density_counts.tif']
  
  
#First run the statistics over voxelized data for volumetric comparison          
  g_control=gs.read_group(control)
  g_fam=gs.read_group(fam)
  g_unfam=gs.read_group(unfam)   


  def t_test(group1, g1, group2, g2, label):       
      pvals, sign=gs.t_test_voxelization(group1, group2, signed = True, 
                                         remove_nan = True, p_cutoff = 0.05)
      
      pvc=gs.color_p_values(pvals, sign, positive = [1,0], negative = [0,1], 
                            p_cutoff = None, positive_trend = None, 
                            negative_trend = None, pmax = None)
      
      g1avg=gs.mean(g1)
      g2avg=gs.mean(g2)
      
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/VoxelStatisticsNeg'+label+'.tif',pvc[:,:,:,0]) #Group2<Group1
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/VoxelStatisticsPos'+label+'.tif',pvc[:,:,:,1]) #Group2>Group1
      
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/AVG_Control'+label+'.tif',g1avg)
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/AVG_Exposed'+label+'.tif',g2avg)

  t_test(group1=control, g1=g_control, group2=fam, g2=g_fam, label='control_vs_fam')
  t_test(group1=control, g1=g_control, group2=unfam, g2=g_unfam, label='control_vs_unfam')
  t_test(group1=fam, g1=g_fam, group2=unfam, g2=g_unfam, label='fam_vs_unfam')
  
 #%%
#Run statistics over brain areas to find the ones with a significant difference
 
  control = ['/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F1Control/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F10Control/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F11Control/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F16Control/cells_ClearMap1_points_transformed.npy']          
            
  fam = ['/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F3Fam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F5Fam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F7Fam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F9Fam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F13Fam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F15Fam/cells_ClearMap1_points_transformed.npy']
  
  unfam = ['/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F2Unfam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F4Unfam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F6Unfam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F8Unfam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F12Unfam/cells_ClearMap1_points_transformed.npy',
            '/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F14Unfam/cells_ClearMap1_points_transformed.npy']

  points_control=gs.read_group(control, combine=False)
  points_fam=gs.read_group(fam, combine=False)
  points_unfam=gs.read_group(unfam, combine=False)
  
  def test_brain_areas(points1, points2):  
      #Set the proper annotation file
      annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(1,256)), orientation=(1,2,3),
      overwrite=False, verbose=True)
  
      counts1=gs.count_points_group_in_regions(points1, annotation_file = annotation_file, 
                                           weight_group = None, invalid = 0, hierarchical = False)
      counts2=gs.count_points_group_in_regions(points2, annotation_file = annotation_file, 
                                           weight_group = None, invalid = 0, hierarchical = False)  
      labels=ano.annotation.names
  
      pvalreg, signreg=gs.t_test_region_countss(counts1, counts2, 
                                                annotation_file = annotation_file, 
                                                signed = True, remove_nan = True,
                                                p_cutoff = 0.05, equal_var = False)
      
  test_brain_areas(points1=points_control, points2=points_fam)  
  test_brain_areas(points1=points_unfam, points2=points_fam)  
  test_brain_areas(points1=points_control, points2=points_unfam)  

  
  
        
