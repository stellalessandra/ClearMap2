# %% Initialize workspace
import sys
sys.path.append('/data/user/ClearMap2')
sys.path.append('/data01/astella/ClearMap2')
import scipy.io
import numpy as np
import numpy.lib.recfunctions as rfn
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import analyze_cells_energy as ace
import matplotlib.pyplot as plt

##############################################################################
# Initialization of parameters
root_directory = '/data01/astella/Projects/'
experiment = 'SexualImprinting'
experimental_group = 'Black_wf_WholeBrain'
data_directory = root_directory + experiment + '/' \
                + experimental_group + '/'
subjects = [name for name in os.listdir(data_directory) \
            if os.path.isdir(os.path.join(data_directory, name))]
# load query file where we added volumes for each area
volumes = pd.read_csv('volumes.csv')
thresholds = [4500, 3500, 2500, 1500, 500]

def load_subject(subject, data_directory, threshold):
    file_suffix = 'cells_' + str(threshold)
    df_mouse = pd.read_csv(data_directory + subject + '/' + file_suffix + '.csv')
    df_mouse = ace.reformat_df_mouse(df=df_mouse)
    
    return df_mouse

for threshold in thresholds:
    print(threshold)
    dict_results_across_mice = {subject: ace.calculate_cells_energy_per_level(
        df_mouse=load_subject(subject=subject,
                            data_directory=data_directory,
                            threshold=threshold),
                            vol=volumes,
                            level=8) for subject in subjects}

    np.save('dict_results/dict_results_across_mice_Black_wf_WholeBrain'+str(threshold)+'.npy', dict_results_across_mice)


