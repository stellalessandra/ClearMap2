#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import analyze_cells_energy as ace
import re
import utils
import itertools
import seaborn as sns

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
    for mouse in dict_results.keys():
        dictionary_labels = {area:volumes.loc[volumes['safe_name'] == area]['acronym'].values[0]\
                             for area in dict_results[mouse]['area'].to_list()}
        dict_results[mouse]['acronym'] = dictionary_labels.values()
    return dict_results

def select_significant_areas(df, threshold=0.05):
    tests = ['pval_Control_vs_Fam', 'pval_Control_vs_Unfam', 'pval_Fam_vs_Unfam']
    d = []
    for key in tests:
        d.append(df.sort_values(by=key)[[
        'area', key]]['area'][df.sort_values(by=key)[[
        'area', key]][key] < threshold].to_list())
    # return flattened list
    return list(itertools.chain.from_iterable(d))

def identify_pls_sig_areas(saliences, threshold, volumes):
    # print areas over the saliences threshold
    threshold1 = threshold
    out1 = saliences.gt(threshold1).apply(lambda x: x.index[x].tolist(), axis=1).to_list()
    threshold2 = -threshold
    out2 = saliences.lt(threshold2).apply(lambda x: x.index[x].tolist(), axis=1).to_list()
    # get area names and not acronyms
    pls_sig_areas = [volumes.loc[volumes['acronym'] == area]['safe_name'].values[0] \
                     for j in range(len(out1)) for area in out1[j]] + [
        volumes.loc[volumes['acronym'] == area]['safe_name'].values[0] \
    for j in range(len(out2)) for area in out2[j]]
    pls_sig_areas = np.unique(np.array(pls_sig_areas))
    return pls_sig_areas


def format_data_pls(dict_results, batch, table):
    data = pd.DataFrame()
    for mouse in dict_results.keys():
        if table=='n_cells':
            temp = dict_results[mouse].filter(['acronym','n_cells'], axis=1).set_index('acronym').T
        elif table=='energy':
            temp = dict_results[mouse].filter(['acronym','energy'], axis=1).set_index('acronym').T
        else:
            raise ValueError('either n_cells or energy input')
        temp.insert(loc=0, column='subject', value=mouse)
        temp.insert(loc=1, column='sex', value='F')
        temp.insert(loc=2, column='group', value=re.split(r'(\d+)', mouse)[-1])
        temp.reset_index(drop = True, inplace = True)
        data = pd.concat([data, temp], axis=0) 
    data = data.loc[:, (data != 0).any(axis=0)]
    if table=='n_cells':
        data.to_csv('./results_pls/'+batch+'_n_cells.csv')
    else:
        data.to_csv('./results_pls/'+batch+'_energy.csv')
        

def create_df_levels(volumes):
    areas_level8 = list(zip(volumes.loc[volumes['st_level'] == 8]['acronym'].values,
                volumes.loc[volumes['st_level'] == 8]['id'].values))
    # remove not analyzed areas at level 8 (areas with lowercase letters)
    areas_level8 = [area for area in areas_level8 if \
                    volumes[volumes['id']==area[1]]['safe_name'].values[0][0].isupper()]

    areas_level5 = list(zip(volumes.loc[volumes['st_level'] == 5]['acronym'].values,
                    volumes.loc[volumes['st_level'] == 5]['id'].values))
    set_ids_level5 = set([area[1] for area in areas_level5])
    # create dataframe with area and corresponding parent at level 5
    df_levels = pd.DataFrame()
    df_levels['area'] = [area[0] for area in areas_level8]
    acronyms_level5 = []
    names_level5 = []
    for area in areas_level8:
        # find id of parent area at st_level 5
        # find set of parents
        # loop over depths and make a set of id
        parents_ids = [int(volumes[volumes['acronym']==area[0]][l].values) \
                         for l in range(1,11) \
                           if not np.isnan(volumes[volumes['acronym']==area[0]][l].values)]
        # identify a parent with st_level 5
        id_parent_level5 = None
        for parent in parents_ids:
            if volumes[volumes['id']==parent]['st_level'].values[0] == 5:
                #identify the parent
                id_parent_level5 = parent
        # find corresponding acronym
        acronyms_level5.append(volumes[volumes['id']==id_parent_level5]['acronym'].values[0])
        # find corresponding name
        names_level5.append(volumes[volumes['id']==id_parent_level5]['safe_name'].values[0])
    df_levels['name_area'] = [volumes[volumes['acronym']==area]['safe_name'].values[0] for area in df_levels['area']]
    df_levels['parent_l5'] = acronyms_level5
    df_levels['name_parent_l5'] = names_level5
#     df_levels.to_csv('areas_parents_level5.csv')
    return df_levels
        

def plot_contrasts(df_data, index, ax):
    df = df_data.iloc[index]
    df = pd.DataFrame(df).T
    df = pd.melt(df)
    sns.barplot(y = df.value, 
        x = df.variable, 
        ax=ax,
        data=df)
    
        
def plot_saliences(df_data, index, ax, df_levels):
    df = df_data.iloc[index]
    df = pd.DataFrame(df).T
    df = pd.melt(df)
    df['Brain Hierarchy'] = [df_levels[df_levels['area']==area]['name_parent_l5'].values[0] for area in df['variable']]
    sns.barplot(y = df.value, 
                x = df.variable, 
                ax=ax, 
                hue=df['Brain Hierarchy'],
                data=df,
                dodge=False)
    
 