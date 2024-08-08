# %% Initialize workspace
import sys
sys.path.append('/data/user/ClearMap2')
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
source = config['source']
size = config['size']
method = config['method']
radius = tuple(config['radius'])
rerun = config['rerun_alignment']
debug_detection = config['debug_detection']
align_to = config['align_to']

# Create test data and debug
slicing = (slice_x, slice_y, slice_z)
thresholds_filt = {
    'source': source,
    'size': size
}

# parsing subject and threshold
parser = argparse.ArgumentParser()
parser.add_argument('subject', 
                    metavar='subject', 
                    type=str,
                    help='mouse which is analyzed in batch file')

parser.add_argument('shape_detection_threshold', 
                    metavar='shape_detection_threshold', 
                    type=int,
                    help='threshold which is analyzed in batch file')
args = parser.parse_args()
subject = args.subject
shape_detection_threshold = args.shape_detection_threshold
shape_detection_threshold = int(args.shape_detection_threshold)

print('subject and threshold', subject, shape_detection_threshold)


if user == 'szucca/Ilaria':
    sys.path.append('/data/szucca/Ilaria')
    data_directory = '/data/szucca/Ilaria/Projects/' + experiment + '/' \
            + experimental_group + '/'+ subject + '/'
else:
    sys.path.append('/data01/' + user)
    data_directory = '/data01/' + user + '/' + experiment + '/' \
                + experimental_group + '/'+ subject + '/'
            
            
# make directories needed for project
if not os.path.exists(data_directory):
    os.makedirs(data_directory)
for folder in ['elastix_auto_to_reference', 
               'elastix_resampled_to_auto', 
               'elastix_resampled_to_reference']:
    if not os.path.exists(data_directory+folder):
        os.makedirs(data_directory+folder)
        
        
# Workspace initialization
#expression_raw = 'cFos/Z<Z,4> .tif'
expression_raw = 'cFos/stitched_fused_tp_0_ch_0.tif'
expression_auto = 'Auto/Z<Z,4> .tif'

ws = wsp.Workspace('CellMap', directory=data_directory)
ws.update(raw=expression_raw, autofluorescence=expression_auto)
ws.info()
ws.debug = False

###########################################################################


# creation resources directory
resources_directory = settings.resources_path


# Cell detection
cell_detection(ws=ws, 
               slicing=slicing, 
               shape=shape_param, 
               threshold_detection=shape_detection_threshold, 
               debugging=debug_detection)

# Cell filtering
cell_filtering(ws=ws, 
               thresholds_filtering=thresholds_filt, 
               threshold_detection=shape_detection_threshold,
               debugging=debug_detection)


# Alignment and annotation of detected and filtered results
cell_alignment_and_annotation(ws=ws, 
                              threshold_detection=shape_detection_threshold,
                              orientation=orientation, 
                              slicing=slicing,
                              align_to=align_to)


# exports of results
export_csv(ws=ws, threshold_detection=shape_detection_threshold)

