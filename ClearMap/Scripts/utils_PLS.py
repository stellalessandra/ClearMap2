#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import analyze_cells_energy as ace
import re
import utils
import itertools
import seaborn as sns
import utils

def heatmap(file, limit, scale, mask, xsize, ysize, zsize, normalize):
    import numpy as np
    import pandas as pd

    axis = np.arange(-limit, limit + 1)
    x, y, z = np.meshgrid(axis, axis, axis)
    ball = np.heaviside((limit ** 2 - (x - limit) ** 2 - (y - limit) ** 2 - (z - limit) ** 2), 1)
    output = np.zeros((zsize, ysize, xsize))
    temp = pd.read_csv(file)
    for index, row in temp.iterrows():
        cx = np.rint(row['x'] * scale).astype(int)
        cy = np.rint(row['y'] * scale).astype(int)
        cz = np.rint(row['z'] * scale).astype(int)
        mx = cx - limit
        my = cy - limit
        mz = cz - limit
        Mx = cx + limit + 1
        My = cy + limit + 1
        Mz = cz + limit + 1
        if (mx >= 0) & (my >= 0) & (mz >= 0) & (Mx < xsize) & (My < ysize) & (Mz < zsize):
            if mask[cz, cy, cx] == 1:
                output[mz:Mz, my:My, mx:Mx] += ball

    if normalize:
        output = (output * 50000) / temp.shape[0]
    output = output * mask

    return output


def read_heat_folder(folder, xsize, ysize, zsize):
    import numpy as np
    from os import listdir
    from os.path import join
    import tifffile as tiff

    lista = listdir(folder)
    array = np.zeros((len(lista), zsize, ysize, xsize))
    n = 0
    for name in lista:
        file = join(folder, name)
        array[n, ...] = tiff.imread(file)
        n = n+1

    return array


def ttest(mean1, mean2, std1, std2, n1, n2):
    import numpy as np
    from scipy import stats

    t = (mean1 - mean2)/np.sqrt((std1**2)/n1 + (std2**2)/n2)
    df = ((std1**2)/n1 + (std2**2)/n2)**2/(((std1**2)/n1)**2/(n1-1) +((std2**2)/n2)**2/(n2-1))

    return 2 - 2*stats.t.cdf(t, df)


def pls(x, y):
    import numpy as np
    diag = np.diagflat((np.ones((1, y.shape[0])) @ y) ** (-1))
    m = diag @ y.transpose() @ x
    r = m - np.ones((m.shape[0], 1)) @ ((np.ones((1, m.shape[0])) @ m) / m.shape[0])
    u, s, v = np.linalg.svd(r, full_matrices=False)
    return u, s, v


def procrustes(u, s, v, u0):
    import numpy as np
    n, o, p = np.linalg.svd(np.matmul(u0.transpose(), u), full_matrices=False)
    q = n @ p.transpose()
    vr = v.transpose() @ q
    ur = u @ np.diagflat(s) @ q
    return ur, vr.transpose()


def bootstrap_test(x, y, v0, u0, n, proc):
    import numpy as np
    vdist = np.zeros((n,) + v0.shape)
    m = x.shape[0]
    for i in np.arange(n):
        # generate random index sequence for bootstrapping (i.e. sampling with replacement)
        while True:
            idx = np.random.randint(0, m, m)
            # extract resampled arrays
            xsh = x[idx]
            ysh = y[idx]
            if not np.any(np.all(ysh[..., :] == 0, axis=0)):
                break
        u, s, v = pls(xsh, ysh)
        if proc:
            ur, vr = procrustes(u, s, v, u0)
        else:
            vr = v
        vdist[i, ...] = vr
    vs = np.std(vdist, axis=0)
    return vs

################################################ utils functions defined by AS

def reformat_dict_acronym(dict_results, volumes):
    """
    Replaces area names in a dictionary with their corresponding acronyms.

    Args:
        dict_results (dict): A dictionary containing results.
        volumes (pandas.DataFrame): DataFrame with volume data.

    Returns:
        dict: The modified dictionary with acronyms.
    """
    for mouse in dict_results.keys():
        dictionary_labels = {area: volumes.loc[volumes['safe_name'] == area]['acronym'].values[0]
                             for area in dict_results[mouse]['area'].to_list()}
        dict_results[mouse]['acronym'] = list(dictionary_labels.values())
    return dict_results


def select_significant_areas(df, threshold=0.05):
    """
    Selects significant areas based on p-values in a DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame.
        threshold (float, optional): Threshold for significance. Defaults to 0.05.

    Returns:
        list: List of significant areas.
    """
    tests = ['pval_Control_vs_Fam', 'pval_Control_vs_Unfam', 'pval_Fam_vs_Unfam']
    d = []
    for key in tests:
        d.append(df.sort_values(by=key)[['area', key]]['area'][df.sort_values(by=key)[['area', key]][key] < threshold].to_list())
    # Return flattened list
    return list(itertools.chain.from_iterable(d))


def identify_pls_sig_areas(saliences, threshold, volumes):
    """
    Identifies areas with significant saliences based on a threshold.

    Args:
        saliences (pandas.DataFrame): DataFrame with salience values.
        threshold (float): Threshold for significance.
        volumes (pandas.DataFrame): DataFrame with volume data.

    Returns:
        numpy.ndarray: Array of significant area names.
    """
    # Print areas over the saliences threshold
    threshold1 = threshold
    out1 = saliences.gt(threshold1).apply(lambda x: x.index[x].tolist(), axis=1).to_list()
    threshold2 = -threshold
    out2 = saliences.lt(threshold2).apply(lambda x: x.index[x].tolist(), axis=1).to_list()

    # Get area names and not acronyms
    pls_sig_areas = [volumes.loc[volumes['acronym'] == area]['safe_name'].values[0]
                     for j in range(len(out1)) for area in out1[j]] + [
        volumes.loc[volumes['acronym'] == area]['safe_name'].values[0]
        for j in range(len(out2)) for area in out2[j]]
    pls_sig_areas = np.unique(np.array(pls_sig_areas))
    return pls_sig_areas



def format_data_pls(dict_results, batch, table):
    """
    Formats data from a dictionary of results into a pandas DataFrame.

    Args:
        dict_results (dict): A dictionary containing mouse-specific data.
        batch (str): The batch identifier.
        table (str): The type of data to extract (e.g., 'n_cells', 'energy', etc.).

    Returns:
        pd.DataFrame: A DataFrame containing formatted data.

    Raises:
        ValueError: If an invalid table name is provided.

    Example:
        dict_results = {
            'mouse1': pd.DataFrame(...),
            'mouse2': pd.DataFrame(...)
        }
        formatted_data = format_data_pls(dict_results, 'BatchA', 'n_cells')
    """
    data = pd.DataFrame()  # Initialize an empty DataFrame to store the formatted data

    # Validate the input table name
    if table not in ['n_cells', 'energy', 'density', 'relative_density']:
        raise ValueError('Invalid input value for "table". Please choose from: n_cells, energy, density, relative_density.')

    # Process data for each mouse
    for mouse in dict_results.keys():
        # Extract relevant columns and transpose the DataFrame
        temp = dict_results[mouse].filter(['acronym', table], axis=1).set_index('acronym').T

        # Add additional columns: subject, sex, and group
        temp.insert(loc=0, column='subject', value=mouse)
        temp.insert(loc=1, column='sex', value='F')  # Assuming all mice are female
        temp.insert(loc=2, column='group', value=utils._split_string(mouse)[-1])  # Extract group information

        # Reset the index and concatenate with existing data
        temp.reset_index(drop=True, inplace=True)
        data = pd.concat([data, temp], axis=0)

    # Remove columns with all zeros
    data = data.loc[:, (data != 0).any(axis=0)]

    return data

        
def create_df_levels(volumes):
    """
    Creates a DataFrame containing information about brain areas and their corresponding parent areas.

    Args:
        volumes (pd.DataFrame): A DataFrame containing brain area data.

    Returns:
        pd.DataFrame: A DataFrame with columns for area, parent at level 5, and related names.

    Example:
        volumes = pd.DataFrame(...)  # Load brain area data
        df_levels = create_df_levels(volumes)
    """
    areas_level8 = list(zip(volumes.loc[volumes['st_level'] == 8]['acronym'].values,
                            volumes.loc[volumes['st_level'] == 8]['id'].values))
    
    # Remove areas with lowercase letters at level 8
    areas_level8 = [area for area in areas_level8 if volumes[volumes['id'] == area[1]]['safe_name'].values[0][0].isupper()]

    areas_level5 = list(zip(volumes.loc[volumes['st_level'] == 5]['acronym'].values,
                            volumes.loc[volumes['st_level'] == 5]['id'].values))
    set_ids_level5 = set([area[1] for area in areas_level5])

    # Create a DataFrame to store area information
    df_levels = pd.DataFrame()
    df_levels['area'] = [area[0] for area in areas_level8]
    acronyms_level5 = []
    names_level5 = []

    for area in areas_level8:
        # Find the parent area ID at st_level 5
        # find set of parents
        # loop over depths and make a set of id
        parents_ids = [int(volumes[volumes['acronym'] == area[0]][l].values) \
                       for l in range(1, 11) \
                       if not np.isnan(volumes[volumes['acronym'] == area[0]][l].values)]
        
        # Identify a parent with st_level 5
        id_parent_level5 = None
        for parent in parents_ids:
            if volumes[volumes['id'] == parent]['st_level'].values[0] == 5:
                id_parent_level5 = parent
        
        # Find corresponding acronym and name
        acronyms_level5.append(volumes[volumes['id'] == id_parent_level5]['acronym'].values[0])
        names_level5.append(volumes[volumes['id'] == id_parent_level5]['safe_name'].values[0])

    df_levels['name_area'] = [volumes[volumes['acronym'] == area]['safe_name'].values[0] for area in df_levels['area']]
    df_levels['parent_l5'] = acronyms_level5
    df_levels['name_parent_l5'] = names_level5

    # Uncomment the following line to save the DataFrame to a CSV file
    # df_levels.to_csv('area_levels.csv', index=False)

    return df_levels

        

def plot_contrasts(df_data, index, ax, palette):
    """
    Plots bar charts for contrasts from a DataFrame of contrasts output of PLS analysis.

    Args:
        df_data (pd.DataFrame): A DataFrame containing contrast data.
        index (int): Index of the specific contrast to plot.
        ax (matplotlib.axes.Axes): The target axes for plotting.
        palette (list or str): Color palette for the bars.

    Returns:
        None

    Example:
        df_data = pd.DataFrame(...)  # Load contrast data
        fig, ax = plt.subplots()
        plot_contrasts(df_data, 0, ax, 'Set1')
    """
    df = df_data.iloc[index]
    df = pd.DataFrame(df).T
    df = pd.melt(df)
    sns.barplot(y=df.value,
                x=df.variable,
                ax=ax,
                data=df,
                palette=palette)

    
        
def plot_saliences(df_data, index, ax, df_levels, palette):
    """
    Plots bar charts for salience values output of PLS analysis with additional hierarchy information.

    Args:
        df_data (pd.DataFrame): A DataFrame containing salience data.
        index (int): Index of the specific salience to plot.
        ax (matplotlib.axes.Axes): The target axes for plotting.
        df_levels (pd.DataFrame): A DataFrame with brain area hierarchy information.
        palette (list or str): Color palette for the bars.

    Returns:
        None

    Example:
        df_data = pd.DataFrame(...)  # Load salience data
        df_levels = pd.DataFrame(...)  # Load brain area hierarchy information
        fig, ax = plt.subplots()
        plot_saliences(df_data, 0, ax, df_levels, 'husl')
    """
    df = df_data.iloc[index]
    df = pd.DataFrame(df).T
    df = pd.melt(df)
    df['Brain Hierarchy'] = [df_levels[df_levels['area'] == area]['name_parent_l5'].values[0] for area in df['variable']]
    sns.barplot(y=df.value,
                x=df.variable,
                ax=ax,
                hue=df['Brain Hierarchy'],
                data=df,
                dodge=False,
                palette=palette)

    

def plot_panel_contrasts(batch, variable, palette='tab10'):
    """
    Plots panel contrasts for a given batch and variable.

    Args:
        batch (str): The batch identifier (e.g., 'batch1', 'batch2').
        variable (str): The variable of interest (e.g., 'score', 'response_time').
        palette (str, optional): Color palette for plotting. Defaults to 'tab10'.

    Returns:
        matplotlib.figure.Figure, numpy.ndarray: The figure containing the contrast plots and an array of axes.
    """
    # Read the contrasts data from a CSV file
    contrasts = pd.read_csv(f'./results_pls/{batch}_{variable}_contrasts.csv')
    
    # Rename columns for better readability
    contrasts = contrasts.rename(columns={col: col[6:] for col in contrasts.columns if col.startswith('group_')})
    
    # number of contrasts
    n_contrasts = contrasts.shape[0]
    
    # Create a 1x3 grid of subplots (one for each contrast)
    fig, axes = plt.subplots(1, n_contrasts, sharey='row', figsize=(10, 5))
    
    # Plot each contrast
    for i in range(n_contrasts):
        plot_contrasts(df_data=contrasts, index=i, ax=axes[i], palette=palette)
        axes[i].tick_params(axis='x', labelrotation=90)
        if i == 0:
            axes[i].set_ylabel('Contrast')  # Set y-label for the first subplot
        else:
            axes[i].set_ylabel('')  # Hide y-label for subsequent subplots
        axes[i].set(xlabel=None)  # Remove x-label
        
    plt.tight_layout()  # Adjust subplot spacing
    return fig, axes

    
    
def plot_panel_saliences(batch, variable, df_levels, palette=sns.color_palette("Paired")):
    """
    Plots panel saliences for a given batch and variable.

    Args:
        batch (str): The batch identifier (e.g., 'batch1', 'batch2').
        variable (str): The variable of interest (e.g., 'score', 'response_time').
        df_levels (pandas.DataFrame): DataFrame containing levels information.
        palette (list, optional): Color palette for plotting. Defaults to seaborn's "Paired".

    Returns:
        matplotlib.figure.Figure, numpy.ndarray: The figure containing the salience plots and an array of axes.
    """
    # Read the saliences data from a CSV file
    saliences = pd.read_csv(f'./results_pls/{batch}_{variable}_saliences.csv')
    
    # number of saliences
    n_saliences = saliences.shape[0]  
    
    # Create a nx1 grid of subplots (one for each salience)
    fig, axes = plt.subplots(n_saliences, 1, sharex='row', figsize=(13, 7))
    plt.subplots_adjust(top=0.9, left=0.03, right=0.8)
     
    
    # Keeps every nth label on the x-axis
    n = 4
    
    # Plot each salience
    for i in range(n_saliences):
        plot_saliences(df_data=saliences, index=i, ax=axes[i], df_levels=df_levels, palette=palette)
        axes[i].tick_params(axis='x', labelrotation=90)
        axes[i].axhline(y=2.57, linestyle='-.', color='darkgrey')
        axes[i].axhline(y=-2.57, linestyle='-.', color='darkgrey')
        
        # Hide some x-axis labels
        [l.set_visible(False) for (i, l) in enumerate(axes[i].xaxis.get_ticklabels()) if i % n != 0]
        
        axes[i].set_ylabel('z-score(salience)')
        
        if i != n_saliences-1:
            axes[i].set_xlabel('')
            axes[i].set(xticklabels=[])
    
    # Set x-label for the last subplot
    axes[n_saliences-1].set_xlabel('Area')
    
    # Add legend to the first subplot
    axes[0].legend(loc='right', bbox_to_anchor=(1.25, 0.3))
    
    # Remove legends from the other two subplots
    for i in range(1,n_saliences):
        axes[i].get_legend().remove()
    
    return fig, axes

    
    
    
    
    