import yaml
from yaml import Loader
import os

def list_subjects():
    try:
        with open("ClearMap/Scripts/configfile.yaml", 'r') as stream:
            config = yaml.load(stream, Loader=Loader)
    except FileNotFoundError:
        with open("configfile.yaml", 'r') as stream:
            config = yaml.load(stream, Loader=Loader) 
    user = config['user']
    experiment = config['experiment']
    experimental_group = config['experimental_group']
    
    data_directory = '/data01/' + user + '/Projects/' + experiment + '/' \
                + experimental_group + '/'
                
    subjects = [name for name in os.listdir(data_directory) \
                if os.path.isdir(os.path.join(data_directory, name))]
    
    return subjects

def divide_in_exp_groups(list_subjects):
    experimental_groups = {'Control':[],
                          'Unfam':[],
                          'Fam':[]}
    for subject in list_subjects:
        for exp_group in experimental_groups.keys():
            if exp_group in subject:
                experimental_groups[exp_group].append(subject)
    return experimental_groups
     