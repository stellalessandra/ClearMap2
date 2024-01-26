import yaml
from yaml import Loader
import os
import re

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


def split_string(input_string):
    # Define the regular expression pattern
    pattern = re.compile(r'([A-Za-z]+[0-9]+)(.*)')

    # Use the pattern to match the input string
    match = pattern.match(input_string)

    # Check if there is a match
    if match:
        # Extract the matched groups
        prefix = match.group(1)
        rest = match.group(2)

        # Return the result as a list
        return [prefix, rest]
    else:
        # If there's no match, return the original string
        return [input_string]

    
def divide_in_exp_groups(list_subjects, group_labels = ['Control', 'Fam', 'Unfam']):
    experimental_groups = {label:[] for label in group_labels}
    for subject in list_subjects:
        for exp_group in experimental_groups.keys():
            if exp_group in split_string(input_string=subject):
                experimental_groups[exp_group].append(subject)
    return experimental_groups
     