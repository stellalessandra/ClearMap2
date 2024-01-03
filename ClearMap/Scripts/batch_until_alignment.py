import os
import yaml
from yaml import Loader
# loading parameters
with open("ClearMap/Scripts/configfile.yaml", 'r') as stream:
    config = yaml.load(stream, Loader=Loader)

subjects = config['subjects']

for subject in subjects:
#    os.system('python ClearMap/Scripts/analysis_after_stitched.py '+ subject)
    os.system('python ClearMap/Scripts/analysis_until_alignment.py '+ subject)
