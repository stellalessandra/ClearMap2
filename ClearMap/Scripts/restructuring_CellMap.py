#%% Initialize workspace
import sys
sys.path.append('/data01/astella/ClearMap2')
from ClearMap.Environment import *  #analysis:ignore
import argparse
import os
import numpy.lib.recfunctions as rfn
import scipy.io

def convert_data_to_numpy(rerun=False):
    if rerun and not os.path.exists(directory+'stitched.npy'):
        print("Converting data to stitched file...")
        source = ws.source('raw');
        # TODO: check if stitched file exists already
        sink   = ws.filename('stitched')
        io.delete_file(sink)
        io.convert(source, sink, processes=None, verbose=False);
        io.mhd.write_header_from_source(ws.filename('stitched'))
    else:
        print("Stitching has already been done!")
    
    
def resampling(source_res, sink_res, directory, rerun=False):
    """
    Parameters
    ----------
    source_res: triplet of x,y,z for source resolution
    sink_res: triplet of x,y,z for sink resolution
    """
    print("Doing resampling step...")
    resample_parameter = {
      "source_resolution" : (source_res[0], source_res[1], source_res[2]),
      "sink_resolution"   : (sink_res[0], sink_res[1], sink_res[2]),
      "processes" : 4,
      "verbose" : False,             
      };
    if rerun and not os.path.exists(directory+'resampled.npy'):    
        io.delete_file(ws.filename('resampled'))
        res.resample(ws.filename('raw'), sink=ws.filename('resampled'), 
                     **resample_parameter)
    else:
        print("Resampling has already been done!")

def initialization_alignment(directory):
    annotation_file, reference_file, distance_file = ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(1,256)), orientation=(1,2,3),
      overwrite=False, verbose=False);
          
    #alignment parameter files    
    align_channels_affine_file   = io.join(directory, 
                                           'Alignment/align_affine.txt')
    align_reference_affine_file  = io.join(directory, 
                                           'Alignment/align_affine.txt')
#    align_reference_bspline_file = io.join(resources_directory, 
#                                           'Alignment/align_bspline.txt')
    align_reference_bspline_file = io.join(directory, 
                                           'Alignment/align_bsplineTest.txt')
    
    return annotation_file, reference_file, distance_file, \
align_channels_affine_file, align_reference_affine_file, \
align_reference_bspline_file
    

def alignment_resampled_to_autofluorescence(align_channels_affine_file, 
                                            directory,
                                            rerun):
    print("Aligning resampled data to autofluorescence file")
    
    if rerun and not os.path.exists(directory+
                                    'elastix_resampled_to_auto/results0.zraw'):
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
    else:
        print("Already done!")
    
    
def alignment_autofluorescence_to_reference(reference_file,
                                            align_reference_affine_file,
                                            align_reference_bspline_file,
                                            directory,
                                            rerun):
    print("Aligning autofluorescence file to reference file")
    if rerun and not os.path.exists(directory+
                                    'elastix_auto_to_reference/results0.zraw'):
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
    else:
        print("Already done!")
    
    
def alignment(rerun=False):
    """
    Script performing the alignment first on resampled data to autofluorescence
    and then to reference file
    """
    print("Doing the alignment...")
    annotation_file, reference_file, distance_file, \
align_channels_affine_file, align_reference_affine_file, \
align_reference_bspline_file = initialization_alignment(resources_directory)

    alignment_resampled_to_autofluorescence(align_channels_affine_file, 
                                            directory,
                                            rerun)
    
    alignment_autofluorescence_to_reference(reference_file,
                                            align_reference_affine_file,
                                            align_reference_bspline_file,
                                            directory,
                                            rerun)


def create_test_data(slicing):    
    #select sublice for testing the pipeline
    ws.create_debug('stitched', slicing=slicing);
    ws.debug = False; 
    io.mhd.write_header_from_source(ws.filename('stitched', prefix='debug'))
    

def visualization_filtering(threshold_detection):
    # visualization of filtering for test debug data
    coordinates = np.hstack([ws.source('cells', 
                                       postfix='filtered_'+str(threshold_detection), 
                                       prefix='debug')[c][:,None] for c in 'xyz']);
    source=ws.source('stitched', prefix='debug')
    aa=np.zeros(source.shape, dtype='int16')
    aa[(coordinates[:,0],coordinates[:,1],coordinates[:,2])]=1;
    io.write(ws.filename('cells', postfix='overlap_debug'), 
             np.asarray(aa, order='F'))
    io.mhd.write_header_from_source(ws.filename('cells', 
                                                postfix='overlap_debug'))
    plt.show()
    
    
def cell_detection_filtering(slicing, shape, threshold_detection,
                             thresholds_filtering, debugging=True):
    cell_detection_parameter = cells.default_cell_detection_parameter.copy();
    cell_detection_parameter['illumination_correction'] = None;
    cell_detection_parameter['background_correction']['shape'] = shape;
#    cell_detection_parameter['background_correction']['save'] = ws.filename('stitched', postfix='background')
    cell_detection_parameter['background_correction']['save'] = None
#    cell_detection_parameter['intensity_detection']['method'] = ['mean'];
    cell_detection_parameter['intensity_detection']['measure'] = ['source'];
    cell_detection_parameter['shape_detection']['threshold'] = threshold_detection;
  
#    cell_detection_parameter['dog_filter']['shape']=None;
  
    io.delete_file(ws.filename('cells', postfix='maxima'))
#    cell_detection_parameter['maxima_detection']['save'] = ws.filename('cells', postfix='maxima')
    cell_detection_parameter['maxima_detection']['save'] = None

#    cell_detection_parameter['maxima_detection']['h_max']=None;

    processing_parameter = cells.default_cell_detection_processing_parameter.copy();
    processing_parameter.update(
        processes = 6, # 'serial',6
        size_max = 20, #100, #35,20
        size_min = 11,# 30, #30,11
        overlap  = 10, #32, #10,
        verbose = True
        )
    
    # doing cell detection
    if debugging == True:
        # creating test data
        create_test_data(slicing)
        # doing cell detection
        cells.detect_cells(ws.filename('stitched', prefix='debug'), 
                           ws.filename('cells', postfix=str(threshold_detection), 
                                       prefix='debug'),
                           cell_detection_parameter=cell_detection_parameter, 
                           processing_parameter=processing_parameter)
        # cell filtering
        cells.filter_cells(source = ws.filename('cells',
                                            postfix=str(threshold_detection), 
                                            prefix='debug'), 
                           sink = ws.filename('cells', 
                                              postfix='filtered'+str(threshold_detection),
                                              prefix='debug'),
                                              thresholds=thresholds_filtering)
        print("plotting...")
        visualization_filtering(threshold_detection)
    else:
        # doing cell detection
        cells.detect_cells(ws.filename('stitched'), 
                             ws.filename('cells', postfix=str(threshold_detection)),
                             cell_detection_parameter=cell_detection_parameter, 
                             processing_parameter=processing_parameter)
        # cell filtering
        cells.filter_cells(source = ws.filename('cells',
                                            postfix=str(threshold_detection)), 
                       sink = ws.filename('cells', 
                                          postfix='filtered'+str(threshold_detection)),
        thresholds=thresholds_filtering)


def visualization_cell_statistics(threshold_detection, directory,
                                  show=True):
      source = ws.source('cells', postfix='filtered'+str(threshold_detection))
      plt.figure(1)
      plt.clf()
      names = source.dtype.names
      nx,ny = p3d.subplot_tiling(len(names))
      for i, name in enumerate(names):
          plt.subplot(nx, ny, i+1)
          plt.hist(source[name])
          plt.title(name)
      plt.tight_layout()
      if show:
          plt.show()
      else:
          if os.path.exists(directory+'figures/'):
              plt.savefig(directory+'figures/filtering_stats.png')
          else: 
              mkdirp(directory+'figures/')
              plt.savefig(directory+'figures/filtering_stats.png')
      # Decide where to save plot
      
      
def transformation(coordinates):
    coordinates = res.resample_points(
                    coordinates, sink=None, orientation=None, 
                    source_shape=io.shape(ws.filename('stitched')), 
                    sink_shape=io.shape(ws.filename('resampled')))
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory=ws.filename('resampled_to_auto'), 
                    binary=True, indices=False)
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory=ws.filename('auto_to_reference'),
                    binary=True, indices=False)
    return coordinates
     

def cell_alignment_and_annotation(threshold_detection, orientation):
      # cell alignment
      source = ws.source('cells', postfix='filtered'+str(threshold_detection))
      coordinates = np.array([source[c] for c in 'xyz']).T
      coordinates_transformed = transformation(coordinates)
      # cell annotation
      # set the common annotation file
      annotation_file, reference_file, \
      distance_file=ano.prepare_annotation_files(
              slicing=(slice(None),slice(None),slice(1,256)), 
              orientation=orientation,
              overwrite=False, verbose=True)
      
      label = ano.label_points(coordinates_transformed,
                               annotation_file=annotation_file, key='order')
      names = ano.convert_label(label, key='order', value='name')
      # save results
      coordinates_transformed.dtype=[(t,float) for t in ('xt','yt','zt')]
      label = np.array(label, dtype=[('order', int)]);
      names = np.array(names, dtype=[('name' , 'U256')])
      cells_data = rfn.merge_arrays([source[:], coordinates_transformed, 
                                     label, names], flatten=True, 
      usemask=False)
  
      io.write(ws.filename('cells'), cells_data)


def export_csv():
    source = ws.source('cells')
    header = ', '.join([h[0] for h in source.dtype.names])
    np.savetxt(ws.filename('cells', extension='csv'), source[:], header=header, 
               delimiter=',', fmt='%s')
    
    
def export_clearmap1():
    print("Exporting ClearMap1 file...")
    source = ws.source('cells')
    clearmap1_format = {'points' : ['x', 'y', 'z'], 
                       'points_transformed' : ['xt', 'yt', 'zt'],
                      ' intensities' : ['source', 'dog', 'background', 'size']}
    for filename, names in clearmap1_format.items():
        sink = ws.filename('cells', postfix=['ClearMap1', filename])
        data = np.array([source[name] if name in source.dtype.names else np.full(source.shape[0], np.nan) for name in names])
        io.write(sink, data)
    
    
def export_matlab(threshold_detection):
    print("Exporting matlab file...")
    matStructure = {}
    values = np.load(ws.filename('cells'))  #load the cells.npy file containing all the saved data
    variable = os.path.splitext('cells.npy')[0]  #remove the .npy extenstion
    variable = variable.lstrip('0123456789.-_ ')  # -- Might not be necessary -- check that all names starts with a number/letter
    matStructure[variable] = values  #load the values in the structure
    filename = directory + 'cells'+str(threshold_detection)+'Mat' + '.mat' #define the location and name of the mat file
    scipy.io.savemat(filename, matStructure) #generate and save the .mat file
    
    
def voxelization(orientation, method='sphere', radius=(7,7,7)):
    print("Voxelization...")
    source = ws.source('cells')
    coordinates = np.array([source[n] for n in ['xt','yt','zt']]).T
    intensities = source['source']
  
    #  To compare data from different animals we need to move coordinates back to the same reference atlas
    annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(1,256)), orientation=orientation,
      overwrite=False, verbose=True)
    
    # Unweighted
    voxelization_parameter = dict(
        shape = io.shape(annotation_file), 
        dtype = None, 
        weights = None,
        method = method, 
        radius = radius, 
        kernel = None, 
        processes = None, 
        verbose = True
      )
    vox.voxelize(coordinates, sink=ws.filename('density', postfix='counts'), 
                 **voxelization_parameter)
    
    # Weighted
    voxelization_parameter = dict(
        shape = io.shape(annotation_file),
        dtype = None, 
        weights = intensities,
        method = method, 
        radius = radius, 
        kernel = None, 
        processes = None, 
        verbose = True
      )
    vox.voxelize(coordinates, 
                 sink=ws.filename('density', postfix='intensities'), 
                 **voxelization_parameter)
  
  

if __name__ == "__main__":
    # parsing name of the folder
    parser = argparse.ArgumentParser(
        description='user')
    parser.add_argument('user',
                        metavar='user',
                        type=str,
                        help='usually in form of initial+surname (username in SuperBrain)')
    args = parser.parse_args()
    user = args.user
    sys.path.append('/data01/'+user)
    directory = '/data01/'+user+'/Projects/SexualImprinting/C57_MaleUrine_Exposure_cFos/F1Control/' 
    
    # TODO: what about this?
    expression_raw      = 'cFos/Z<Z,4> .tif'         
    expression_auto     = 'Auto/Z<Z,4> .tif'
      
    ws = wsp.Workspace('CellMap', directory=directory);
    ws.update(raw=expression_raw, autofluorescence=expression_auto)
    ws.info()
    ws.debug = False
      
    # creation resources directory
    resources_directory = settings.resources_path
    
    # convertion of data to numpy
    convert_data_to_numpy()
    
    # resampling of autofluorescence
    resampling(source_res=[1.626,1.626,5], sink_res=[25,25,25], 
               directory=resources_directory)
    
    # alignment of resampled to autofluorescence and to reference
    alignment(rerun=False)
    
    # Create test data and debug
    slicing = (slice(2900,3100),slice(4120,4320),slice(620,650))
    shape_param = (7,7)
    shape_detection_threshold = 700
    thresholds_filt = {
      'source' : (10,5000),
      'size'   : (20,100)
      }
    orientation = (1,2,3)
    
    cell_detection_filtering(slicing=slicing, shape=shape_param, 
                             threshold_detection=shape_detection_threshold,
                             thresholds_filtering=thresholds_filt, 
                             debugging=True)
    
    visualization_cell_statistics(threshold_detection=shape_detection_threshold,
                                  directory=resources_directory)
    
    cell_alignment_and_annotation(threshold_detection=shape_detection_threshold, 
                                  orientation=orientation)
    
    # exports
    export_csv()
    export_clearmap1()
    export_matlab(threshold_detection=shape_detection_threshold)
    
    # voxelization
    voxelization(orientation, method='sphere', radius=(7,7,7))
  
    
    
    
    
    
    
    
