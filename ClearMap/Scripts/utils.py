import yaml
from yaml import Loader
import os
import re


def _split_string(input_string):
    """
    This function splits a string into two parts based on a regular expression pattern.
    
    Args:
        input_string (str): The string to be split.
        
    Returns:
        list: A list containing the two parts of the split string.
    """
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
    

def list_subjects():
    """
    This function reads a configuration file and returns a list of subjects 
    from a specific directory.
    
    Returns:
        subjects (list): A list of subject names.
    """
    # Try to open the configuration file in the specified directory
    try:
        with open("ClearMap/Scripts/configfile_ila.yaml", 'r') as stream:
            config = yaml.load(stream, Loader=Loader)
    # If the file is not found in the specified directory, look for it in the current directory
    except FileNotFoundError:
        with open("configfile_ila.yaml", 'r') as stream:
            config = yaml.load(stream, Loader=Loader) 
    
    # Extract the data directory path from the configfile
    data_directory = config['data_directory']
                
    # Get a list of all directories (subjects) in the data directory
    subjects = [name for name in os.listdir(data_directory) \
                if os.path.isdir(os.path.join(data_directory, name))]
    
    return subjects

    
def divide_in_exp_groups(list_subjects, group_labels = ['Control', 'Fam', 'Unfam']):
    """
    This function divides subjects into experimental groups based on their names.
    
    Args:
        list_subjects (list): A list of subject names.
        group_labels (list, optional): A list of experimental group labels. 
            Defaults to ['Control', 'Fam', 'Unfam'].
            
    Returns:
        dict: A dictionary where the keys are the experimental group labels and 
            the values are lists of subjects belonging to each group.
    """
    # Initialize a dictionary to hold the experimental groups
    experimental_groups = {label:[] for label in group_labels}
    
    # Loop over the subjects
    for subject in list_subjects:
        # Loop over the experimental groups
        for exp_group in experimental_groups.keys():
            # If the experimental group label is in the subject name, add the subject to the group
            if exp_group in _split_string(input_string=subject):
                experimental_groups[exp_group].append(subject)
    
    return experimental_groups

    
def split_density_maps(path):
    """
    Split an image in half vertically and process it for the right hemisphere.

    This function loads a TIFF image, cuts it in half vertically keeping the second half,
    swaps the x and z axes, inverts the x-axis, and saves the resulting image.

    Parameters:
    path (str): The directory path where the 'density_counts.tif' image is located and where
                the processed image 'density_counts_right_hemisphere.tif' will be saved.
    """
    # Load the image from the specified path
    image = imread(path + 'density_counts.tif')

    # Define the index to split the image vertically and keep the second half
    half_index = 201
    image = image[:, :, half_index:]

    # Swap the x and z axes of the image
    image = np.swapaxes(image, 0, 2)

    # Invert the x-axis
    image = image[::-1, :, :]

    # Save the processed image to the specified path
    imsave(path + 'density_counts_right_hemisphere.tif', image)
     
