# %% Initialize workspace
import sys
sys.path.append('/data01/astella/ClearMap2')
from ClearMap.Environment import *  # analysis:ignore
import scipy.io
import numpy as np
import numpy.lib.recfunctions as rfn
import os
import argparse
import yaml
from yaml import Loader
import time
from ClearMap.Scripts.functions_CellMap import *

##############################################################################
# Initialization of parameters


print('Starting the analysis', time.time())
initial_time = time.time()
times = []

with open("ClearMap/Scripts/configfile.yaml", 'r') as stream:
    config = yaml.load(stream, Loader=Loader)

user = config['user']
experiment = config['experiment']
experimental_group = config['experimental_group']
source_res = config['source_res']
sink_res = config['sink_res']
orientation = tuple(config['orientation'])
slice_x = config['slice_x']
slice_y = config['slice_y']
slice_z = config['slice_z']
shape_param = config['shape_param']
shape_detection_threshold = config['shape_detection_threshold']
source = config['source']
size = config['size']
method = config['method']
radius = tuple(config['radius'])
rerun = config['rerun']
debug_detection = config['debug_detection']

# Create test data and debug
slicing = (slice_x, slice_y, slice_z)
thresholds_filt = {
    'source': source,
    'size': size
}

# parsing name of the folder
parser = argparse.ArgumentParser(
    description='subject')
parser.add_argument(
    'subject',
    metavar='subject',
    type=str,
    help='mouse which is analyzed in batch file')
args = parser.parse_args()
subject = args.subject


sys.path.append('/data01/' + user)
# makedir here with subject
directory = '/data01/' + user + '/Projects/' + experiment + '/' \
            + experimental_group + '/'+ subject + '/'
            
            
# make directories needed for project
if not os.path.exists(directory):
    os.makedirs(directory)
for folder in ['Auto', 'cFos', 'elastix_auto_to_reference', 
               'elastix_resampled_to_auto']:
    if not os.path.exists(directory+folder):
        os.makedirs(directory+folder)
        
        
# Workspace initialization
expression_raw = 'cFos/Z<Z,4> .tif'
expression_auto = 'Auto/Z<Z,4> .tif'

ws = wsp.Workspace('CellMap', directory=directory)
ws.update(raw=expression_raw, autofluorescence=expression_auto)
ws.info()
ws.debug = False

###########################################################################


# creation resources directory
resources_directory = settings.resources_path



# convertion of data to numpy
convert_data_to_numpy(ws=ws, directory=directory, rerun=rerun)



# resampling of autofluorescence
resampling(ws=ws, source_res=source_res, sink_res=sink_res,
           directory=resources_directory)

print('Resampling done', time.time() - initial_time)
times.append(time.time() - initial_time)



# alignment of resampled to autofluorescence and to reference
alignment(ws=ws, rerun=rerun)

print('Alignment done', time.time() - times[-1])
times.append(time.time() - times[-1])



# Cell detection and filtering
cell_detection_filtering(ws=ws, slicing=slicing, shape=shape_param,
                         threshold_detection=shape_detection_threshold,
                         thresholds_filtering=thresholds_filt,
                         debugging=debug_detection)

print('Detection and filtering done', time.time() - times[-1])
times.append(time.time() - times[-1])



# Visualization cell statistics
visualization_cell_statistics(ws=ws, 
    threshold_detection=shape_detection_threshold,
    directory=resources_directory)



# Alignment and annotation of detected and filtered results
cell_alignment_and_annotation(ws=ws, 
    threshold_detection=shape_detection_threshold,
    orientation=orientation)

print('Cell alignment and annotation done', time.time() - times[-1])
times.append(time.time() - times[-1])



# exports of results
export_csv(ws=ws)
export_clearmap1(ws=ws)
export_matlab(ws=ws, threshold_detection=shape_detection_threshold,
              directory=directory)



# voxelization
voxelization(ws=ws, orientation=orientation, method=method, radius=radius)

print('Detection and filtering done', time.time() - times[-1])
times.append(time.time() - times[-1])
np.save('ClearMap/Scripts/times'+subject, times)
np.save(directory+'params', config)
