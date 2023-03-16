import yaml
from yaml import Loader
import os

def list_subjects():
    with open("ClearMap/Scripts/configfile.yaml", 'r') as stream:
        config = yaml.load(stream, Loader=Loader)
    
    user = config['user']
    experiment = config['experiment']
    experimental_group = config['experimental_group']
    
    data_directory = '/data01/' + user + '/Projects/' + experiment + '/' \
                + experimental_group + '/'
                
    subjects = [name for name in os.listdir(data_directory) \
                if os.path.isdir(os.path.join(data_directory, name))]
    
    return subjects