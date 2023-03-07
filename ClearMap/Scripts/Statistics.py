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
  
  def t_test(group1, group2, label1, label2): 
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
      
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/VoxelStatisticsNeg-'+label1+'_'+label2+'.tif',pvc[:,:,:,0]) #Group2<Group1
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/VoxelStatisticsPos-'+label1+'_'+label2+'.tif',pvc[:,:,:,1]) #Group2>Group1
      
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/AVG_'+label1+'.tif',g1avg)
      io.write('/data01/astella/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/AVG_'+label2+'.tif',g2avg)

#  t_test(group1=control, group2=fam, label1='control', label2='fam')
#  t_test(group1=control, group2=unfam, label1='control', label2='unfam')
#  t_test(group1=fam, group2=unfam, label1='fam', label2='unfam')
  
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
      return pvalreg, signreg, labels
      
  pvalreg_control_fam, signreg_control_fam, labels = test_brain_areas(points1=points_control, points2=points_fam)  
  pvalre_unfam_fam, signreg_unfam_fam, labels = test_brain_areas(points1=points_unfam, points2=points_fam)  
  pvalreg_control_unfam, signreg_control_unfam, labels = test_brain_areas(points1=points_control, points2=points_unfam)  

  
  
        
