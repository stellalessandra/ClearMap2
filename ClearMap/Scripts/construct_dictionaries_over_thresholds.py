# !pip install odfpy
# !pip install seaborn
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import yaml
from yaml import Loader
import numpy as np
import analyze_cells_energy as ace
import utils
from scipy.stats import ttest_ind, mannwhitneyu
import networkx as nx
import matplotlib.pyplot as plt

root_directory = '/home/stella/Documents/Torino/projects/'
experiment = 'SexualImprinting'
experimental_group = 'Black_wf_WholeBrain'
data_directory = root_directory + experiment + '/' \
                + experimental_group + '/'
subjects = [name for name in os.listdir(data_directory) \
            if os.path.isdir(os.path.join(data_directory, name))]

# load query file where we added volumes for each area
volumes = ace.clean_volumes_database()

def load_subject(subject, data_directory, threshold):
    file_suffix = 'cells_' + str(threshold)
    df_mouse = pd.read_csv(data_directory + subject + '/' + file_suffix + '.csv')
    df_mouse = ace.reformat_df_mouse(df=df_mouse)
    
    return df_mouse

thresholds = [500, 1500, 2500, 3500, 4500]

for subject in subjects:
    for threshold in thresholds:
            result = ace.calculate_cells_energy_per_level(
                 df_mouse=load_subject(subject=subject,
                                       data_directory=data_directory,
                                       threshold=threshold),
                                       vol=volumes,
                                       level=8)
            np.save(data_directory+subject+'/'+subject+str(threshold)+'.npy', result)