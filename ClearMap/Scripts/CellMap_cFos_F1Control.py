#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CellMap
=======

This script is the main pipeline to analyze immediate early gene expression 
data from iDISCO+ cleared tissue [Renier2016]_.

See the :ref:`CellMap tutorial </CellMap.ipynb>` for a tutorial and usage.


.. image:: ../Static/cell_abstract_2016.jpg
   :target: https://doi.org/10.1016/j.cell.2020.01.028
   :width: 300

.. figure:: ../Static/CellMap_pipeline.png

  iDISCO+ and ClearMap: A Pipeline for Cell Detection, Registration, and 
  Mapping in Intact Samples Using Light Sheet Microscopy.


References
----------
.. [Renier2016] `Mapping of brain activity by automated volume analysis of immediate early genes. Renier* N, Adams* EL, Kirst* C, Wu* Z, et al. Cell. 2016 165(7):1789-802 <https://doi.org/10.1016/j.cell.2016.05.007>`_
"""
__author__    = 'Christoph Kirst <christoph.kirst.ck@gmail.com>'
__license__   = 'GPLv3 - GNU General Pulic License v3 (see LICENSE)'
__copyright__ = 'Copyright © 2020 by Christoph Kirst'
__webpage__   = 'http://idisco.info'
__download__  = 'http://www.github.com/ChristophKirst/ClearMap2'

if __name__ == "__main__":
     
  #%%############################################################################
  ### Initialization 
  ###############################################################################
  
  #%% Initialize workspace
  import sys
  sys.path.append('/data01/astella')
  
  from ClearMap.Environment import *  #analysis:ignore
  
  #directories and files
  directory = '/data01/astella/Projects/SexualImprinting/Black_wf_WholeBrain/BL45_Control/'    
  
  expression_raw      = 'cFos/Z<Z,4> .tif'         
  expression_auto     = 'Auto/Z<Z,4> .tif'
  
  ws = wsp.Workspace('CellMap', directory=directory);
  ws.update(raw=expression_raw, autofluorescence=expression_auto)
  ws.info()
  
  ws.debug = False
  
  resources_directory = settings.resources_path
  

  

  #%%############################################################################
  ### Data conversion
  ############################################################################### 
  
  #%% Convet raw data to npy file     
               
  source = ws.source('raw');
  sink   = ws.filename('stitched')
  io.delete_file(sink)
  io.convert(source, sink, processes=None, verbose=True);
  
  io.mhd.write_header_from_source(ws.filename('stitched'))
  
  
  #%%############################################################################
  ### Resampling and atlas alignment 
  ###############################################################################
        
  #%% Resample 
             
  resample_parameter = {
      "source_resolution" : (1.626, 1.626, 5),
      "sink_resolution"   : (25,25,25),
      "processes" : 4,
      "verbose" : True,             
      };
  
  io.delete_file(ws.filename('resampled'))
  
  res.resample(ws.filename('raw'), sink=ws.filename('resampled'), **resample_parameter)
  
  
  #%% Resample autofluorescence
      
  resample_parameter_auto = {
      "source_resolution" : (1.626, 1.626, 10),
      "sink_resolution"   : (25,25,25),
      "processes" : 4,
      "verbose" : True,                
      };    
 # io.delete_file(ws.filename('resampled', postfix='autofluorescence'))
  
  res.resample(ws.filename('autofluorescence'), sink=ws.filename('resampled', postfix='autofluorescence'), **resample_parameter_auto)
  
  #%% Initialize alignment 
  
  annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(None)), orientation=(3,2,1),
      overwrite=False, verbose=True);
          
  #alignment parameter files    
  align_channels_affine_file   = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_affine_file  = io.join(resources_directory, 'Alignment/align_affine.txt')
#  align_reference_bspline_file = io.join(resources_directory, 'Alignment/align_bspline.txt')
  align_reference_bspline_file = io.join(resources_directory, 'Alignment/align_bsplineTest.txt')
  
  #%% Aignment - resampled to autofluorescence
  
  # align the two channels
  align_channels_parameter = {            
      #moving and reference images
      "moving_image" : ws.filename('resampled', postfix='autofluorescence'),
#      "moving_image" : ws.filename('resampled'),
      "fixed_image"  : ws.filename('resampled'),
      
      #elastix parameter files for alignment
      "affine_parameter_file"  : align_channels_affine_file,
      "bspline_parameter_file" : None,
      
      #directory of the alig'/home/nicolas.renier/Documents/ClearMap_Ressources/Par0000affine.txt',nment result
      "result_directory" :  ws.filename('resampled_to_auto')
      }; 
  
  elx.align(**align_channels_parameter);
  
  #%% Alignment - autoflourescence to reference
  
  # align autofluorescence to reference
  align_reference_parameter = {            
      #moving and reference images
      "moving_image" : reference_file,
      "fixed_image"  : ws.filename('resampled', postfix='autofluorescence'),
#      "fixed_image"  : ws.filename('resampled'),
      
      #elastix parameter files for alignment
      "affine_parameter_file"  :  align_reference_affine_file,
      "bspline_parameter_file" :  align_reference_bspline_file,
      #directory of the alignment result
      "result_directory" :  ws.filename('auto_to_reference')
      };
  
  elx.align(**align_reference_parameter);
  
  
  #%%############################################################################
  ### Create test data
  ###############################################################################
  
  #%% Crop test data 
  
  #select sublice for testing the pipeline
  slicing = (slice(2900,3100),slice(4120,4320),slice(620,650));
  ws.create_debug('stitched', slicing=slicing);
  ws.debug = False; 
  io.mhd.write_header_from_source(ws.filename('stitched', prefix='debug'))
    
  
  #%%############################################################################
  ### Cell detection
  ###############################################################################
  
  #%% Cell detection for debugging:
  
  cell_detection_parameter = cells.default_cell_detection_parameter.copy();
  cell_detection_parameter['illumination_correction'] = None;
  cell_detection_parameter['background_correction']['shape'] = (7,7);
#  cell_detection_parameter['background_correction']['save'] = ws.filename('stitched', postfix='background')
  cell_detection_parameter['background_correction']['save'] = None
#  cell_detection_parameter['intensity_detection']['method'] = ['mean'];
  cell_detection_parameter['intensity_detection']['measure'] = ['source'];
  cell_detection_parameter['shape_detection']['threshold'] = 700;
  
#  cell_detection_parameter['dog_filter']['shape']=None;
  
  io.delete_file(ws.filename('cells', postfix='maxima'))
#  cell_detection_parameter['maxima_detection']['save'] = ws.filename('cells', postfix='maxima')
  cell_detection_parameter['maxima_detection']['save'] = None

#  cell_detection_parameter['maxima_detection']['h_max']=None;
  
  
  processing_parameter = cells.default_cell_detection_processing_parameter.copy();
  processing_parameter.update(
      processes = 6, # 'serial',6
      size_max = 20, #100, #35,20
      size_min = 11,# 30, #30,11
      overlap  = 10, #32, #10,
      verbose = True
      )
  
  cells.detect_cells(ws.filename('stitched', prefix='debug'), ws.filename('cells', postfix='700', prefix='debug'),
                     cell_detection_parameter=cell_detection_parameter, 
                     processing_parameter=processing_parameter)
  
#  io.mhd.write_header_from_source(ws.filename('cells', postfix='maxima'))
#  io.mhd.write_header_from_source(ws.filename('stitched', postfix='background'))
#  io.mhd.write_header_from_source(ws.filename('stitched', postfix='illumination'))
  
   #%% Filter cells
  
  thresholds = {
      'source' : (10,5000),
      'size'   : (20,100)
      }
  
  cells.filter_cells(source = ws.filename('cells',postfix='700', prefix='debug'), 
                     sink = ws.filename('cells', postfix='filtered700', prefix='debug'), 
                     thresholds=thresholds);

  #%% visualization
  
#  #%%
  # TODO: understand this!
  #coordinates = np.hstack([ws.source('cells', postfix='filtered_debug')[c][:,None] for c in 'xyz']);
  
  coordinates = np.hstack([ws.source('cells', postfix='filtered_700', prefix='debug')[c][:,None] for c in 'xyz']);
  source=ws.source('stitched', prefix='debug')
  aa=np.zeros(source.shape, dtype='int16')
  aa[(coordinates[:,0],coordinates[:,1],coordinates[:,2])]=1;
  io.write(ws.filename('cells',postfix='overlap_debug'), np.asarray(aa, order='F'))
  io.mhd.write_header_from_source(ws.filename('cells', postfix='overlap_debug'))
  
  
  #%% Cell detection for the entire dataset:
  
  cell_detection_parameter = cells.default_cell_detection_parameter.copy();
  cell_detection_parameter['illumination_correction'] = None;
  cell_detection_parameter['background_correction']['shape'] = (7,7);
#  cell_detection_parameter['background_correction']['save'] = ws.filename('stitched', postfix='background')
  cell_detection_parameter['background_correction']['save'] = None
#  cell_detection_parameter['intensity_detection']['method'] = ['mean'];
  cell_detection_parameter['intensity_detection']['measure'] = ['source'];
  cell_detection_parameter['shape_detection']['threshold'] = 700;
  
#  cell_detection_parameter['dog_filter']['shape']=None;
  
  io.delete_file(ws.filename('cells', postfix='maxima'))
#  cell_detection_parameter['maxima_detection']['save'] = ws.filename('cells', postfix='maxima')
  cell_detection_parameter['maxima_detection']['save'] = None

#  cell_detection_parameter['maxima_detection']['h_max']=None;
  
  
  processing_parameter = cells.default_cell_detection_processing_parameter.copy();
  processing_parameter.update(
      processes = 6, # 'serial',6
      size_max = 20, #100, #35,20
      size_min = 11,# 30, #30,11
      overlap  = 10, #32, #10,
      verbose = True
      )
  
  cells.detect_cells(ws.filename('stitched'), ws.filename('cells', postfix='700'),
                     cell_detection_parameter=cell_detection_parameter, 
                     processing_parameter=processing_parameter)
  
#  io.mhd.write_header_from_source(ws.filename('cells', postfix='maxima'))
#  io.mhd.write_header_from_source(ws.filename('stitched', postfix='background'))
#  io.mhd.write_header_from_source(ws.filename('stitched', postfix='illumination'))
  
   #%% Filter cells
  
  thresholds = {
      'source' : (10,5000),
      'size'   : (20,100)
      }
  
  cells.filter_cells(source = ws.filename('cells',postfix='700'), 
                     sink = ws.filename('cells', postfix='filtered700'), 
                     thresholds=thresholds);

 #%% Cell statistics
  
  source = ws.source('cells', postfix='filtered700')
  
  plt.figure(1); plt.clf();
  names = source.dtype.names;
  nx,ny = p3d.subplot_tiling(len(names));
  for i, name in enumerate(names):
    plt.subplot(nx, ny, i+1)
    plt.hist(source[name]);
    plt.title(name)
  plt.tight_layout();
  
  
  #%%############################################################################
  ### Cell atlas alignment and annotation
  ###############################################################################
  
  #%% Cell alignment
  
  source = ws.source('cells', postfix='filtered700')
  
  def transformation(coordinates):
    coordinates = res.resample_points(
                    coordinates, sink=None, orientation=None, 
                    source_shape=io.shape(ws.filename('stitched')), 
                    sink_shape=io.shape(ws.filename('resampled')));
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory=ws.filename('resampled_to_auto'), 
                    binary=True, indices=False);
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory=ws.filename('auto_to_reference'),
                    binary=True, indices=False);
        
    return coordinates;
    
  
  coordinates = np.array([source[c] for c in 'xyz']).T;
  
  coordinates_transformed = transformation(coordinates);
  
  #We want to have the coordinates for all files in the same reference orientation
  #To do so we need to flip the coordinates depending on how the images have been acquired.
  
#  coordinates_transformed[:,0]=-coordinates_transformed[:,0]+2*160 #THIS NEED TO CHANGE DEPENDING ON THE BRAIN ORIENTATION
#  coordinates_transformed[:,1]=coordinates_transformed[:,1]+60 #THIS NEED TO CHANGE DEPENDING ON THE BRAIN ORIENTATION
  
  #%% Cell annotation
  
  #Set the common annotation file
  annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(1,256)), orientation=(1,2,3),
      overwrite=False, verbose=True);
          
  label = ano.label_points(coordinates_transformed,annotation_file=annotation_file, key='order');
  names = ano.convert_label(label, key='order', value='name');
  
  #%% Save results
  
  coordinates_transformed.dtype=[(t,float) for t in ('xt','yt','zt')]
  label = np.array(label, dtype=[('order', int)]);
  names = np.array(names, dtype=[('name' , 'U256')])
  
  import numpy.lib.recfunctions as rfn
  cells_data = rfn.merge_arrays([source[:], coordinates_transformed, label, names], flatten=True, usemask=False)
  
  io.write(ws.filename('cells'), cells_data)
  
  
  
  #%%############################################################################
  ### Cell csv generation for external analysis
  ###############################################################################
  
  #%% CSV export
  
  source = ws.source('cells');
  header = ', '.join([h[0] for h in source.dtype.names]);
  np.savetxt(ws.filename('cells', extension='csv'), source[:], header=header, delimiter=',', fmt='%s')
  
  #%% ClearMap 1.0 export
  
  source = ws.source('cells');
  
  clearmap1_format = {'points' : ['x', 'y', 'z'], 
                      'points_transformed' : ['xt', 'yt', 'zt'],
                      'intensities' : ['source', 'dog', 'background', 'size']}
  
  for filename, names in clearmap1_format.items():
    sink = ws.filename('cells', postfix=['ClearMap1_700', filename]);
    data = np.array([source[name] if name in source.dtype.names else np.full(source.shape[0], np.nan) for name in names]);
    io.write(sink, data);
  
 #%% Save Mat File -- This is to save a .mat file which can be opened with Matlab as a cell array
 import scipy.io  #Import library
 
 matStructure = {}
 values=np.load(ws.filename('cells'))  #load the cells.npy file containing all the saved data
 variable=os.path.splitext('cells.npy')[0]  #remove the .npy extenstion
 variable = variable.lstrip('0123456789.-_ ')  # -- Might not be necessary -- check that all names starts with a number/letter
 
 matStructure[variable] = values  #load the values in the structure
 filename = directory + 'cells700Mat' + '.mat' #define the location and name of the mat file

 scipy.io.savemat(filename, matStructure) #generate and save the .mat file
  
  #%%############################################################################
  ### Voxelization - cell density
  ###############################################################################
  
  source = ws.source('cells')
  
  coordinates = np.array([source[n] for n in ['xt','yt','zt']]).T;
  intensities = source['source'];
  
#  To compare data from different animals we need to move coordinates back to the same reference atlas
  annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(1,256)), orientation=(1,2,3),
      overwrite=False, verbose=True);  
  
  #%% Unweighted 
  
  voxelization_parameter = dict(
        shape = io.shape(annotation_file), 
        dtype = None, 
        weights = None,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink=ws.filename('density', postfix='counts'), **voxelization_parameter);
  
  #%% Weighted 
  
  voxelization_parameter = dict(
        shape = io.shape(annotation_file),
        dtype = None, 
        weights = intensities,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink=ws.filename('density', postfix='intensities'), **voxelization_parameter);
  
  #%%
  # TOFIX
  p3d.plot(ws.filename('density', postfix='intensities'))
