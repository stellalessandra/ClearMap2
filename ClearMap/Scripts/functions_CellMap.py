# %% Initialize workspace
#import sys
#sys.path.append('/data01/astella/ClearMap2')
from ClearMap.Environment import *  # analysis:ignore
import scipy.io
import numpy as np
import numpy.lib.recfunctions as rfn
import os
import argparse
import yaml
from yaml import Loader
import time

def reformat_slicing_parameter(slicing):
    """
    Reformat a slicing parameter for use in array slicing.

    This function takes a slicing parameter, which is a tuple of length 3, where each element
    can either be a tuple indicating the start and end points for slicing along that dimension,
    or any other type indicating no specific slicing (i.e., take the whole dimension).
    It returns a tuple of slice objects corresponding to the slicing parameters.

    Parameters:
    slicing (tuple): A tuple of length 3, where each element can be:
                     - A tuple of two integers indicating the start and end points for slicing.
                     - Any other type indicating no specific slicing (equivalent to slice(None)).

    Returns:
    tuple: A tuple of slice objects for each dimension.
    """
    
    # Re-format slicing parameter for the x dimension
    if type(slicing[0]) is tuple:
        slice_x = slice(slicing[0][0], slicing[0][1])
    else:
        slice_x = slice(None)

    # Re-format slicing parameter for the y dimension
    if type(slicing[1]) is tuple:
        slice_y = slice(slicing[1][0], slicing[1][1])
    else:
        slice_y = slice(None)

    # Re-format slicing parameter for the z dimension
    if type(slicing[2]) is tuple:
        slice_z = slice(slicing[2][0], slicing[2][1])
    else:
        slice_z = slice(None)

    # Combine the slice objects into a tuple and return
    slicing = (slice_x, slice_y, slice_z)
    return slicing



def convert_data_to_numpy(ws, directory, rerun=False):
    """
    Convert data to a stitched NumPy file if it doesn't already exist or if rerun is specified.

    This function checks if a stitched NumPy file exists in the specified directory.
    If the file does not exist or if rerun is set to True, it converts the raw data
    to a stitched file. Otherwise, it prints a message indicating that the stitching
    has already been done.

    Parameters:
    ws: Workspace object providing methods to access data sources and filenames.
    directory (str): The directory where the stitched file is stored.
    rerun (bool): Whether to rerun the conversion even if the stitched file already exists.

    Returns:
    None
    """
    
    # Check if rerun is True or the stitched file does not exist
    if rerun or not os.path.exists(directory + 'stitched.npy'):
        print("Converting data to stitched file...")
        
        # Get the source of the raw data from the workspace
        source = ws.source('raw')
        
        # Get the filename for the stitched data from the workspace
        sink = ws.filename('stitched')
        
        # Optionally delete the existing stitched file
        # io.delete_file(sink)
        
        # Convert the raw data to the stitched file format
        io.convert(source, sink, processes=None, verbose=False)
        
        # Write the header for the stitched file from the source
        io.mhd.write_header_from_source(ws.filename('stitched'))
    else:
        print("Stitching has already been done!")
        


def resampling(ws, source_res, sink_res, align_to, directory, rerun=False):
    """
    Perform resampling of image data to specified resolutions.

    This function resamples image data from a source resolution to a sink resolution.
    It also handles a special case for 'cfos_auto' alignment, adjusting the z-resolution
    differently. The function will resample the data only if rerun is True or if the 
    resampled file does not already exist.

    Parameters:
    ws: Workspace object providing methods to access filenames.
    source_res (tuple): A triplet (x, y, z) specifying the source resolution.
    sink_res (tuple): A triplet (x, y, z) specifying the sink resolution.
    align_to (str): A string indicating special alignment cases, such as 'cfos_auto'.
    directory (str): The directory where the resampled file is stored.
    rerun (bool): Whether to rerun the resampling even if the resampled file already exists.

    Returns:
    None
    """
    
    print("Doing resampling step...")
    
    # Define resampling parameters
    resample_parameter = {
        "source_resolution": (source_res[0], source_res[1], source_res[2]),
        "sink_resolution": (sink_res[0], sink_res[1], sink_res[2]),
        "processes": 4,
        "verbose": False,
    }
    
    # Define resampling parameters for 'cfos_auto' alignment case
    if align_to == 'cfos_auto':
        resample_parameter_auto = {
            "source_resolution": (source_res[0], source_res[1], source_res[2] * 2),
            "sink_resolution": (sink_res[0], sink_res[1], sink_res[2]),
            "processes": 4,
            "verbose": False,
        }
    
    # Check if rerun is True or the resampled file does not exist
    if rerun or not os.path.exists(directory + 'resampled.tif'):
        # Optionally delete the existing resampled file
        io.delete_file(ws.filename('resampled'))
        
        # Perform the resampling
        res.resample(ws.filename('stitched'), sink=ws.filename('resampled'), **resample_parameter)
        
        # Perform resampling for 'cfos_auto' alignment case
        if align_to == 'cfos_auto':
            res.resample(ws.filename('autofluorescence'), 
                         sink=ws.filename('resampled', 
                                          postfix='autofluorescence'), 
                         **resample_parameter_auto)
    else:
        print("Resampling has already been done!")


def initialization_alignment(alignment_files_directory, orientation, slicing):
    """
    Initialize alignment by preparing annotation files and defining alignment parameter file paths.

    This function reformats the slicing parameter, prepares annotation files, and sets up paths for
    various alignment parameter files needed for the alignment process.

    Parameters:
    alignment_files_directory (str): The directory where alignment parameter files are stored.
    orientation (str): The orientation of the slices (e.g., 'coronal', 'sagittal', etc.).
    slicing (tuple): A tuple of length 3 defining the slicing parameters for x, y, and z dimensions.

    Returns:
    tuple: A tuple containing paths to the annotation file, reference file, distance file, and
           alignment parameter files (affine and bspline).
    """
    
    # Re-format slicing parameter
    slicing = reformat_slicing_parameter(slicing)
    
    # Prepare annotation files with the specified slicing and orientation
    annotation_file, reference_file, distance_file = ano.prepare_annotation_files(
        slicing=slicing,
        orientation=orientation,
        overwrite=False, 
        verbose=False)

    # Define paths for alignment parameter files
    align_channels_affine_file = io.join(alignment_files_directory, 'Alignment/align_affine.txt')
    align_reference_affine_file = io.join(alignment_files_directory, 'Alignment/align_affine.txt')
    align_reference_bspline_file = io.join(alignment_files_directory, 'Alignment/align_bspline.txt')
    # Optionally use a different bspline file
    # align_reference_bspline_file = io.join(directory, 'Alignment/align_bsplineTest.txt')

    return annotation_file, reference_file, distance_file, 
            align_channels_affine_file, align_reference_affine_file, 
            align_reference_bspline_file


def alignment_resampled_to_autofluorescence(ws, align_channels_affine_file, directory, rerun):
    """
    Align resampled data to autofluorescence data using affine transformation.

    This function aligns the resampled data to autofluorescence data using specified
    alignment parameters. The alignment is performed using the elastix tool. The function
    will only perform the alignment if rerun is True or if the result file does not exist.

    Parameters:
    ws: Workspace object providing methods to access filenames.
    align_channels_affine_file (str): Path to the affine parameter file used for alignment.
    directory (str): The directory where the alignment result is stored.
    rerun (bool): Whether to rerun the alignment even if the result file already exists.

    Returns:
    None
    """
    
    print("Aligning resampled data to autofluorescence file")
    
    # Check if rerun is True or the result file does not exist
    if rerun or not os.path.exists(directory + 'elastix_resampled_to_auto/result.0.zraw'):
        
        # Define alignment parameters
        align_channels_parameter = {
            # Moving and reference images
            "moving_image": ws.filename('resampled', postfix='autofluorescence'),
            # Optionally use the 'resampled' image as the moving image
            # "moving_image": ws.filename('resampled'),
            "fixed_image": ws.filename('resampled'),
            
            # Elastix parameter files for alignment
            "affine_parameter_file": align_channels_affine_file,
            "bspline_parameter_file": None,
            
            # Result directory
            "result_directory": ws.filename('resampled_to_auto')
        }
        
        # Perform the alignment using elastix
        elx.align(**align_channels_parameter)
    else:
        print("Already done!")


def alignment_autofluorescence_to_reference(ws, reference_file, 
                                            align_reference_affine_file, 
                                            align_reference_bspline_file, 
                                            directory, rerun):
    """
    Align autofluorescence data to a reference file using affine and bspline transformations.

    This function aligns the autofluorescence data to a reference file using specified alignment
    parameters. The alignment is performed using the elastix tool. The function will only perform 
    the alignment if rerun is True or if the result file does not exist.

    Parameters:
    ws: Workspace object providing methods to access filenames.
    reference_file (str): Path to the reference file used for alignment.
    align_reference_affine_file (str): Path to the affine parameter file used for alignment.
    align_reference_bspline_file (str): Path to the bspline parameter file used for alignment.
    directory (str): The directory where the alignment result is stored.
    rerun (bool): Whether to rerun the alignment even if the result file already exists.

    Returns:
    None
    """
    
    print("Aligning autofluorescence file to reference file")
    
    # Check if rerun is True or the result file does not exist
    if rerun or not os.path.exists(directory + 'elastix_auto_to_reference/result.0.zraw'):
        
        # Define alignment parameters
        align_reference_parameter = {
            # Moving and reference images
            "moving_image": reference_file,
            "fixed_image": ws.filename('resampled', postfix='autofluorescence'),
            # Optionally use the 'resampled' image as the fixed image
            # "fixed_image": ws.filename('resampled'),
            
            # Elastix parameter files for alignment
            "affine_parameter_file": align_reference_affine_file,
            "bspline_parameter_file": align_reference_bspline_file,
            
            # Result directory
            "result_directory": ws.filename('auto_to_reference')
        }
        
        # Perform the alignment using elastix
        elx.align(**align_reference_parameter)
    else:
        print("Already done!")
        
        
def alignment_resampled_to_reference(ws, reference_file, 
                                     align_reference_affine_file, 
                                     align_reference_bspline_file, 
                                     directory, rerun):
    """
    Align resampled data to a reference file using affine and bspline transformations.

    This function aligns the resampled data to a reference file using specified alignment
    parameters. The alignment is performed using the elastix tool. The function will only perform 
    the alignment if rerun is True or if the result file does not exist.

    Parameters:
    ws: Workspace object providing methods to access filenames.
    reference_file (str): Path to the reference file used for alignment.
    align_reference_affine_file (str): Path to the affine parameter file used for alignment.
    align_reference_bspline_file (str): Path to the bspline parameter file used for alignment.
    directory (str): The directory where the alignment result is stored.
    rerun (bool): Whether to rerun the alignment even if the result file already exists.

    Returns:
    None
    """
    
    print("Aligning resampled data to reference file")
    
    # Check if rerun is True or the result file does not exist
    if rerun or not os.path.exists(directory + 'elastix_resampled_to_reference/result.0.zraw'):
        
        # Define alignment parameters
        align_reference_parameter = {
            # Moving and reference images
            "moving_image": reference_file,
            "fixed_image": ws.filename('resampled'),
            # Optionally use the 'resampled' image as the fixed image
            # "fixed_image": ws.filename('resampled'),
            
            # Elastix parameter files for alignment
            "affine_parameter_file": align_reference_affine_file,
            "bspline_parameter_file": align_reference_bspline_file,
            
            # Result directory
            "result_directory": ws.filename('resampled_to_reference')
        }
        
        # Perform the alignment using elastix
        elx.align(**align_reference_parameter)
    else:
        print("Already done!")


def alignment(ws, directory, alignment_files_directory, orientation, align_to, slicing, rerun=False):
    """
    Perform alignment of data first to autofluorescence and then to a reference file.

    This function orchestrates the alignment process:
    - Initializes alignment parameters based on orientation and slicing.
    - Performs alignment either directly to a reference file or via autofluorescence, based on `align_to` parameter.

    Parameters:
    ws: Workspace object providing methods to access filenames.
    directory (str): Directory path where alignment results are stored.
    alignment_files_directory (str): Directory containing alignment parameter files.
    orientation (str): Orientation of the data.
    align_to (str): Alignment target, either 'cfos' or 'cfos_auto'.
    slicing (tuple): Slicing parameters for data alignment.
    rerun (bool, optional): Whether to force rerun of alignment if results exist. Defaults to False.

    Returns:
    None
    Raises:
    ValueError: If `align_to` is not 'cfos' or 'cfos_auto'.
    """
    
    print("Doing the alignment...")
    
    # Initialize alignment parameters
    annotation_file, reference_file, distance_file, \
        align_channels_affine_file, align_reference_affine_file, \
        align_reference_bspline_file = initialization_alignment(
                alignment_files_directory=alignment_files_directory, 
                orientation=orientation, 
                slicing=slicing)
    
    # Perform alignment based on align_to parameter
    if align_to == 'cfos':
        alignment_resampled_to_reference(ws=ws,
                                         reference_file=reference_file,
                                         align_reference_affine_file=align_reference_affine_file,
                                         align_reference_bspline_file=align_reference_bspline_file,
                                         directory=directory,
                                         rerun=rerun)
    elif align_to == 'cfos_auto':
        # Align resampled data to autofluorescence
        alignment_resampled_to_autofluorescence(ws=ws,
                                                align_channels_affine_file=align_channels_affine_file,
                                                directory=directory,
                                                rerun=rerun)

        # Then align autofluorescence data to reference
        alignment_autofluorescence_to_reference(ws=ws, 
                                                reference_file=reference_file,
                                                align_reference_affine_file=align_reference_affine_file,
                                                align_reference_bspline_file=align_reference_bspline_file,
                                                directory=directory,
                                                rerun=rerun)
    else:
        # Raise an error if align_to is neither 'cfos' nor 'cfos_auto'
        raise ValueError('Align either to cfos or cfos_auto')


def create_test_data(ws, slicing):
    """
    Create test data by selecting a subset of the stitched data based on slicing parameters.

    This function reformats the slicing parameter and selects a sub-slice of the stitched data 
    for testing purposes. It then writes the header from the source to the debug file.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    slicing (tuple): Slicing parameters for selecting a subset of the data.

    Returns:
    None
    """
    
    # Re-format slicing parameter using reformat_slicing_parameter function
    slicing = reformat_slicing_parameter(slicing)
    
    # Create debug data by selecting a sub-slice of the stitched data
    ws.create_debug('stitched', slicing=slicing)
    
    # Disable debug mode
    ws.debug = False
    
    # Write header from source to the debug file
    io.mhd.write_header_from_source(ws.filename('stitched', prefix='debug'))



def visualization_filtering(ws, threshold_detection):
    """
    Visualize filtering for test debug data.

    This function visualizes the filtering process for test debug data by creating a visualization
    of filtered cell coordinates overlaid on debug data. It writes the visualization to a file and 
    displays it using matplotlib.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    threshold_detection (float): Threshold value used for filtering.

    Returns:
    None
    """
    
    # Get coordinates of filtered cells from debug data
    coordinates = np.hstack([ws.source('cells',
                                       postfix='filtered' + str(threshold_detection),
                                       prefix='debug')[c][:, None] for c in 'xyz'])
    
    # Get source data for visualization
    source = ws.source('stitched', prefix='debug')
    
    # Create an array to store the visualization
    aa = np.zeros(source.shape, dtype='int16')
    
    # Mark filtered cell coordinates in the visualization array
    aa[(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2])] = 1
    
    # Write visualization array to a file
    io.write(ws.filename('cells', postfix='overlap_debug'),
             np.asarray(aa, order='F'))
    
    # Write header for the visualization file
    io.mhd.write_header_from_source(ws.filename('cells',
                                                postfix='overlap_debug'))
    
    # Display the visualization using matplotlib
    plt.show()


def cell_detection(ws, slicing, shape, threshold_detection, debugging=True):
    """
    Perform cell detection using specified parameters.

    This function configures cell detection parameters, including background correction, intensity detection,
    shape detection, and maxima detection. It then performs cell detection either with or without debugging,
    creating test data if debugging is enabled.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    slicing (tuple): Slicing parameters for selecting a subset of the data for testing.
    shape (tuple): Shape parameters for background correction.
    threshold_detection (float): Threshold value for shape detection.
    debugging (bool, optional): Whether to enable debugging mode. Defaults to True.

    Returns:
    None
    """
    
    # Configure cell detection parameters
    cell_detection_parameter = cells.default_cell_detection_parameter.copy()
    cell_detection_parameter['illumination_correction'] = None
    cell_detection_parameter['background_correction']['shape'] = shape
    cell_detection_parameter['background_correction']['save'] = None
    cell_detection_parameter['intensity_detection']['measure'] = ['source']
    cell_detection_parameter['shape_detection']['threshold'] = threshold_detection
    cell_detection_parameter['maxima_detection']['save'] = None

    # Configure processing parameters
    processing_parameter = cells.default_cell_detection_processing_parameter.copy()
    processing_parameter.update(
        processes=6,
        size_max=20,
        size_min=11,
        overlap=10,
        verbose=False
    )

    # Perform cell detection
    if debugging:
        # Create test data for debugging
        create_test_data(ws=ws, slicing=slicing)
        # Perform cell detection with debug data
        cells.detect_cells(
            ws.filename('stitched', prefix='debug'),
            ws.filename('cells', postfix=str(threshold_detection), prefix='debug'),
            cell_detection_parameter=cell_detection_parameter,
            processing_parameter=processing_parameter
        )
    else:
        print('Doing the detection...')
        # Perform cell detection without debug data
        cells.detect_cells(
            ws.filename('stitched'),
            ws.filename('cells', postfix=str(threshold_detection)),
            cell_detection_parameter=cell_detection_parameter,
            processing_parameter=processing_parameter
        )
            
            
def cell_filtering(ws, thresholds_filtering, threshold_detection, debugging=True):
    """
    Perform cell filtering based on specified thresholds.

    This function filters cells based on detection thresholds either with or without debugging,
    generating visualizations if debugging is enabled.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    thresholds_filtering (dict): Thresholds for filtering cells.
    threshold_detection (float): Threshold value used for cell detection.
    debugging (bool, optional): Whether to enable debugging mode. Defaults to True.

    Returns:
    None
    """
    
    # Perform cell detection
    if debugging:
        # Perform cell filtering with debug data
        cells.filter_cells(
            source=ws.filename('cells', postfix=str(threshold_detection), prefix='debug'),
            sink=ws.filename('cells', postfix='filtered' + str(threshold_detection), prefix='debug'),
            thresholds=thresholds_filtering
        )
        print("plotting...")
        # Visualize filtering results
        visualization_filtering(ws, threshold_detection)
    else:
        print('Doing the filtering...')
        # Perform cell filtering without debug data
        cells.filter_cells(
            source=ws.filename('cells', postfix=str(threshold_detection)),
            sink=ws.filename('cells', postfix='filtered' + str(threshold_detection)),
            thresholds=thresholds_filtering
        )



def visualization_cell_statistics(ws, threshold_detection, directory, show=False):
    """
    Visualize statistics of filtered cells.

    This function plots histograms of various statistics from filtered cell data and optionally
    saves the plot to a file.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    threshold_detection (float): Threshold value used for cell detection.
    directory (str): Directory path to save the visualization.
    show (bool, optional): Whether to display the plot. Defaults to False.

    Returns:
    None
    """
    
    # Get filtered cell data
    source = ws.source('cells', postfix='filtered' + str(threshold_detection))
    
    # Create a new figure
    plt.figure(1)
    plt.clf()
    
    # Get column names from the data
    names = source.dtype.names
    
    # Calculate subplot layout
    nx, ny = p3d.subplot_tiling(len(names))
    
    # Plot histograms for each statistic
    for i, name in enumerate(names):
        plt.subplot(nx, ny, i + 1)
        plt.hist(source[name])
        plt.title(name)
    
    # Adjust subplot layout
    plt.tight_layout()
    
    # Display or save the plot
    if show:
        plt.show()
    else:
        # Ensure directory exists for saving
        if not os.path.exists(directory + 'figures/'):
            os.makedirs(directory + 'figures/')
        
        # Save the plot as an image
        plt.savefig(directory + 'figures/filtering_stats.png')



def transformation(ws, coordinates, align_to):
    """
    Perform transformation of coordinates based on alignment.

    This function transforms coordinates based on the specified alignment scenario (cfos_auto or not).
    It resamples points and applies transformations accordingly.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    coordinates (numpy array): Coordinates to be transformed.
    align_to (str): Alignment scenario ('cfos_auto' or other).

    Returns:
    numpy array: Transformed coordinates.
    """
    
    # Resample points to match the target shape
    coordinates = res.resample_points(
        coordinates,
        sink=None,
        orientation=None,
        source_shape=io.shape(ws.filename('stitched')),
        sink_shape=io.shape(ws.filename('resampled'))
    )

    # Perform transformations based on alignment scenario
    if align_to == 'cfos_auto':
        # Transform coordinates to autofluorescence reference
        coordinates = elx.transform_points(
            coordinates,
            sink=None,
            transform_directory=ws.filename('resampled_to_auto'),
            binary=True,
            indices=False
        )

        # Transform coordinates to final reference
        coordinates = elx.transform_points(
            coordinates,
            sink=None,
            transform_directory=ws.filename('auto_to_reference'),
            binary=True,
            indices=False
        )
    else:
        # Transform coordinates directly to reference
        coordinates = elx.transform_points(
            coordinates,
            sink=None,
            transform_directory=ws.filename('resampled_to_reference'),
            binary=True,
            indices=False
        )
    
    return coordinates


def cell_alignment_and_annotation(ws, threshold_detection, orientation, align_to, slicing):
    """
    Perform cell alignment and annotation.

    This function aligns cell coordinates, annotates them, and saves the annotated results.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    threshold_detection (float): Threshold value used for cell detection.
    orientation (str): Orientation parameter for annotation.
    align_to (str): Alignment scenario ('cfos_auto' or other).
    slicing (tuple): Slicing parameters for selecting a subset of the data.

    Returns:
    None
    """
    
    # Re-format slicing parameter
    print('slicing', slicing)
    slicing = reformat_slicing_parameter(slicing)
    
    # Cell alignment
    source = ws.source('cells', postfix='filtered' + str(threshold_detection))
    coordinates = np.array([source[c] for c in 'xyz']).T  # Extract x, y, z coordinates from source
    coordinates_transformed = transformation(ws, coordinates, align_to)  # Transform coordinates based on alignment
    
    # Cell annotation
    # Prepare annotation files
    annotation_file, reference_file, distance_file = ano.prepare_annotation_files(
        slicing=slicing,
        orientation=orientation,
        overwrite=False,
        verbose=False
    )
    print('coordinates_transformed', coordinates_transformed)
    print('annotation file', annotation_file)
    
    # Label points based on transformed coordinates
    label = ano.label_points(coordinates_transformed, annotation_file=annotation_file, key='order')
    names = ano.convert_label(label, key='order', value='name')
    
    # Save results
    coordinates_transformed.dtype = [(t, float) for t in ('xt', 'yt', 'zt')]  # Define data types for transformed coordinates
    label = np.array(label, dtype=[('order', int)])  # Define data type for labels
    names = np.array(names, dtype=[('name', 'U256')])  # Define data type for names
    cells_data = rfn.merge_arrays([source[:], coordinates_transformed, label, names], flatten=True, usemask=False)  # Merge data arrays
    io.write(ws.filename('cells'), cells_data)  # Write merged data to file



def export_csv(ws, threshold_detection):
    """
    Export cell data to a CSV file.

    This function exports cell data stored in the workspace (`ws`) to a CSV file.
    It retrieves data from the 'cells' source, optionally filtered by a threshold value.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    threshold_detection (float): Threshold value used for filtering cell data.

    Returns:
    None
    """
    
    # Retrieve cell data from workspace
    source = ws.source('cells')
    
    # Create header for CSV file
    header = ', '.join([h[0] for h in source.dtype.names]) + ', l1, l2, l3'
    
    # Save data to CSV file
    np.savetxt(ws.filename('cells', postfix=str(threshold_detection), extension='csv'),  # Filename for CSV
               source[:],  # Data to be saved
               header=header,  # Header for CSV file
               delimiter=',',  # Delimiter used in CSV
               fmt='%s')  # Format for saving data



def export_clearmap1(ws):
    """
    Export cell data to ClearMap1 format files.

    This function exports cell data stored in the workspace (`ws`) to multiple ClearMap1 format files.
    It retrieves data from the 'cells' source and formats it according to specified ClearMap1 file formats.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.

    Returns:
    None
    """
    
    print("Exporting ClearMap1 file...")
    
    # Retrieve cell data from workspace
    source = ws.source('cells')
    
    # Define ClearMap1 format specifications
    clearmap1_format = {
        'points': ['x', 'y', 'z'],
        'points_transformed': ['xt', 'yt', 'zt'],
        'intensities': ['source', 'dog', 'background', 'size']
    }
    
    # Iterate over ClearMap1 format items and export data
    for filename, names in clearmap1_format.items():
        sink = ws.filename('cells', postfix=['ClearMap1', filename])  # Define filename for ClearMap1 file
        data = np.array([source[name] if name in source.dtype.names else np.full(source.shape[0], np.nan)
                         for name in names])  # Extract data according to ClearMap1 format
        io.write(sink, data)  # Write data to ClearMap1 format file


def export_matlab(ws, threshold_detection, directory):
    """
    Export cell data to a MATLAB (.mat) file.

    This function exports cell data stored in the workspace (`ws`) to a MATLAB (.mat) file.
    It loads cell data from a numpy (.npy) file, formats it into a MATLAB structure, and saves it.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    threshold_detection (float): Threshold value used for filtering cell data.
    directory (str): Directory path where the MATLAB file will be saved.

    Returns:
    None
    """
    
    print("Exporting MATLAB file...")
    
    # Load cell data from numpy file
    values = np.load(ws.filename('cells'))
    
    # Remove file extension and strip unwanted characters from variable name
    variable = os.path.splitext('cells.npy')[0].lstrip('0123456789.-_ ')
    
    # Create MATLAB structure
    matStructure = {variable: values}
    
    # Define filename for MATLAB (.mat) file
    filename = directory + 'cells' + str(threshold_detection) + 'Mat' + '.mat'
    
    # Save data to MATLAB file
    scipy.io.savemat(filename, matStructure)



def voxelization(ws, orientation, slicing, method='sphere', radius=(7, 7, 7)):
    """
    Perform voxelization of cell coordinates and intensities.

    This function voxelizes cell coordinates and intensities from the workspace (`ws`).
    It prepares annotation files, sets voxelization parameters, and saves voxelized data
    into density files for counts and intensities.

    Parameters:
    ws: Workspace object providing methods to access filenames and data.
    orientation (str): Orientation of the data.
    slicing (tuple): Slicing parameters for data subsetting.
    method (str, optional): Voxelization method ('sphere' or 'cube'). Defaults to 'sphere'.
    radius (tuple, optional): Radius of the voxelization in each dimension. Defaults to (7, 7, 7).

    Returns:
    None
    """
    
    # Re-format slicing parameter
    slicing = reformat_slicing_parameter(slicing)
    print("Voxelization...")
    
    # Retrieve cell data from workspace
    source = ws.source('cells')
    coordinates = np.array([source[n] for n in ['xt', 'yt', 'zt']]).T
    intensities = source['source']

    # Prepare annotation files for voxelization
    annotation_file, reference_file, distance_file = ano.prepare_annotation_files(
        slicing=slicing, orientation=orientation,
        overwrite=False, verbose=True)

    # Voxelization parameters for unweighted counts
    voxelization_parameter = {
        'shape': io.shape(annotation_file),
        'dtype': None,
        'weights': None,
        'method': method,
        'radius': radius,
        'kernel': None,
        'processes': None,
        'verbose': True
    }
    
    # Perform voxelization for counts
    vox.voxelize(coordinates, sink=ws.filename('density', postfix='counts'),
                 **voxelization_parameter)

    # Voxelization parameters for weighted intensities
    voxelization_parameter['weights'] = intensities
    
    # Perform voxelization for intensities
    vox.voxelize(coordinates,
                 sink=ws.filename('density', postfix='intensities'),
                 **voxelization_parameter)

