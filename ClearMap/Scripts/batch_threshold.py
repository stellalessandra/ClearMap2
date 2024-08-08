import os
import yaml
from yaml import Loader
# loading parameters
with open("ClearMap/Scripts/configfile.yaml", 'r') as stream:
    config = yaml.load(stream, Loader=Loader)

subjects = config['subjects']
thresholds = config['shape_detection_threshold']

for subject in subjects:
    for threshold in thresholds:
        os.system('python ClearMap/Scripts/analysis_threshold.py '+subject +' '+str(threshold))
