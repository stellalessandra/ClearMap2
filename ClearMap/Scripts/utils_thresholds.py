import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import yaml
from yaml import Loader
import numpy as np
import analyze_cells_energy as ace
import utils
import matplotlib.colors as cm
import utils_PLS as upls


def sum_values(subjects, thresholds, data_directory):
    """
    Summarize various metrics (n_cells, energy, density, relative_density) for given subjects and thresholds.

    This function reads CSV files for each subject and threshold, sums the specified metrics, and compiles the
    results into a single DataFrame.

    Parameters:
    subjects (list): List of subjects (e.g., ['subject1', 'subject2']).
    thresholds (list): List of thresholds (e.g., [1, 2, 3]).
    data_directory (str): Directory path where the subject data files are located.

    Returns:
    pd.DataFrame: A DataFrame containing summed values for each metric, subject, and threshold.
    """
    
    # Initialize an empty DataFrame with the specified columns
    data = pd.DataFrame(columns=['subject', 'threshold', 'n_cells', 'energy', 'density', 'relative_density'])
    
    # Create a list of subjects repeated for each threshold
    data['subject'] = [subject for subject in subjects for threshold in thresholds]
    
    # Create a list of thresholds repeated for each subject
    data['threshold'] = [threshold for subject in subjects for threshold in thresholds]
    
    # Iterate through each metric and calculate the sum for each subject and threshold
    for metric in ['n_cells', 'energy', 'density', 'relative_density']:
        data[metric] = [pd.read_csv(f"{data_directory}/{subject}/{subject}{threshold}.csv")[metric].sum()
                        for subject in subjects for threshold in thresholds]
    
    return data


def plot_values(data, palette='Set2'):
    """
    Plot various metrics (n_cells, energy, density, relative_density) for given subjects and thresholds.

    This function creates a 2x2 grid of line plots for the specified metrics. Each line plot shows the values
    of a metric across subjects, colored by threshold.

    Parameters:
    data (pd.DataFrame): DataFrame containing columns 'subject', 'threshold', 'n_cells', 'energy', 'density', and 'relative_density'.
    palette (str): Color palette for the line plots.

    Returns:
    None
    """
    # Create a 2x2 grid of subplots
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    plt.subplots_adjust(hspace=0.06)  # Adjust vertical space between plots

    # Plot for 'n_cells' metric
    sns.lineplot(data=data, x='subject', y='n_cells', hue='threshold', ax=axes[0][0], 
                 legend=False, palette=palette, marker="o")
    axes[0][0].tick_params(labelbottom=False)  # Remove x-axis labels for this plot
    axes[0][0].set_xlabel('')  # Remove x-axis label

    # Plot for 'energy' metric
    sns.lineplot(data=data, x='subject', y='energy', hue='threshold', ax=axes[0][1], 
                 palette=palette, marker="o")
    axes[0][1].legend().set_title('')  # Remove legend title
    axes[0][1].tick_params(labelbottom=False)  # Remove x-axis labels for this plot
    axes[0][1].set_xlabel('')  # Remove x-axis label

    # Plot for 'density' metric
    sns.lineplot(data=data, x='subject', y='density', hue='threshold', ax=axes[1][0], 
                 legend=False, palette=palette, marker="o")
    axes[1][0].tick_params(axis='x', rotation=90)  # Rotate x-axis labels for readability
    axes[1][0].set_xlabel('')  # Remove x-axis label

    # Plot for 'relative_density' metric
    sns.lineplot(data=data, x='subject', y='relative_density', hue='threshold', ax=axes[1][1], 
                 legend=False, palette=palette, marker="o")
    axes[1][1].tick_params(axis='x', rotation=90)  # Rotate x-axis labels for readability
    axes[1][1].set_xlabel('')  # Remove x-axis label

    # Display the plot
    plt.show()

    
def find_significant_areas(dictionary, experimental_groups, groups=['Control', 'Fam', 'Unfam'], 
                           test='mannwhitneyu'):
    """
    Find significant brain areas based on non-parametric statistical tests.

    This function calculates cell counts across experimental groups, performs a non-parametric test 
    (Mann-Whitney U by default) to compare groups, and identifies significant brain areas by 
    evaluating p-values. Significant areas are determined based on a p-value threshold (< 0.05).

    Parameters:
    dictionary (dict): A dictionary containing cell count data for each subject.
    experimental_groups (dict): A dictionary mapping experimental groups to their respective subjects.
    groups (list of str): A list of group names to compare. Default is ['Control', 'Fam', 'Unfam'].
    test (str): The statistical test to use. Default is 'mannwhitneyu'.

    Returns:
    pd.DataFrame: A DataFrame listing significant brain areas for each comparison.
    """
    # Calculate cell counts across experimental groups
    dfs = ace.calculate_value_across_groups(experimental_groups=experimental_groups, 
                                            dict_results_across_mice=dictionary, 
                                            value='n_cells')
    
    # Perform statistical tests across groups
    df_ttest = ace.test_across_groups(dfs, groups=groups, test=test)
    
    # Identify significant areas based on p-values
    columns = df_ttest.loc[:, df_ttest.columns != 'area'].columns
    df_sigareas = pd.DataFrame()
    for col in columns:
        df = df_ttest.sort_values(by=col)[['area', col]]
        # Select areas with p-value < 0.05
        significant_areas = df[df[col] < 0.05]['area'].reset_index(drop=True)
        df_sigareas = pd.concat([df_sigareas, significant_areas], axis=1)
    
    # Rename columns to reflect the comparisons
    df_sigareas.columns = [col.replace('pval_', '') for col in columns]
    
    return df_sigareas


def calculate_ratio_detected_cells(subjects, thresholds, data_directory):
    """
    Calculate the ratio of analyzed cells to detected cells for given subjects and thresholds.

    This function computes the ratio of the number of cells analyzed to the number of cells 
    detected for a list of subjects across various thresholds. It uses the summed cell counts 
    from the `sum_values` function and reads the detected cell counts from CSV files.

    Parameters:
    subjects (list): List of subject identifiers.
    thresholds (list): List of threshold values used for analysis.

    Returns:
    pd.DataFrame: A DataFrame containing the subject, threshold, number of cells analyzed,
                  number of detected cells, and the ratio of analyzed to detected cells.
    """
    # Initialize an empty DataFrame to store the results
    data = pd.DataFrame(columns=['subject', 'threshold', 'n_cells_analysis', 'n_detected_cells', 'ratio'])

    # Get the summed values for each subject and threshold
    data_sum = sum_values(subjects=subjects, thresholds=thresholds, data_directory=data_directory)

    # Populate the DataFrame with subject and threshold information
    data['subject'] = [subject for subject in subjects for threshold in thresholds]
    data['threshold'] = [threshold for subject in subjects for threshold in thresholds]

    # Initialize lists to store the analyzed and detected cell counts, and their ratio
    n_cells = []
    n_detected_cells = []
    ratio = []

    # Calculate the number of analyzed and detected cells, and the ratio
    for subject in subjects:
        for threshold in thresholds:
            # Get the total number of analyzed cells for the given subject and threshold
            sum_n = data_sum[(data_sum['subject'] == subject) & 
                             (data_sum['threshold'] == threshold)]['n_cells'].item()
            n_cells.append(sum_n)

            # Read the number of detected cells from the CSV file
            det_n = len(pd.read_csv(data_directory + subject + '/cells_' + str(threshold) + '.csv',
                                    low_memory=False))
            n_detected_cells.append(det_n)

            # Calculate the ratio of analyzed to detected cells
            ratio.append(sum_n / det_n)

    # Store the calculated values in the DataFrame
    data['n_cells_analysis'] = n_cells
    data['n_detected_cells'] = n_detected_cells
    data['ratio'] = ratio

    return data


def dataframe_areas(area, subjects, thresholds, data_directory):
    """
    Create a DataFrame for a specific brain area across different subjects and thresholds.

    This function constructs a DataFrame containing various metrics (n_cells, energy, density, 
    relative_density) for a specified brain area. The data is collected for multiple subjects 
    and thresholds by reading corresponding CSV files.

    Parameters:
    area (str): The specific brain area of interest.
    subjects (list): A list of subject identifiers.
    thresholds (list): A list of threshold values used for analysis.

    Returns:
    pd.DataFrame: A DataFrame with columns 'subject', 'threshold', 'n_cells', 'energy', 
                  'density', and 'relative_density', containing data for the specified area.
    """
    # Initialize the DataFrame with appropriate columns
    data = pd.DataFrame(columns=['subject', 'threshold', 'n_cells', 'energy', 'density', 'relative_density'])

    # Populate the 'subject' and 'threshold' columns
    data['subject'] = [subject for subject in subjects for threshold in thresholds]
    data['threshold'] = [threshold for subject in subjects for threshold in thresholds]

    # Iterate over the metrics and collect data
    for metric in ['n_cells', 'energy', 'density', 'relative_density']:
        values = []
        for subject in subjects:
            for threshold in thresholds:
                # Read the CSV file for the current subject and threshold
                v = pd.read_csv(data_directory + subject + '/' + subject + str(threshold) + '.csv')
                # Extract the value for the specified area and metric
                values.append(v[v['area'] == area][metric].values[0])
        # Assign the collected values to the DataFrame for the current metric
        data[metric] = values

    return data


def corr_matrix(experimental_groups,
               dict_results_across_mice,
               value='n_cells'):
    n_groups = len(experimental_groups.keys())
    dfs = ace.calculate_value_across_groups(experimental_groups=experimental_groups, 
                                  dict_results_across_mice=dict_results_across_mice, 
                                  value=value)
    # rename areas with acronyms
    for df in dfs:
        for i in df.index:
            df.at[i, 'area'] = volumes[volumes['safe_name'] == \
                                       df.at[i, 'area']]['acronym'].values[0]
    # remove null values
    for df in dfs:
        df = df.set_index('area').loc[
                ~(df.set_index('area')==0).all(axis=1)].dropna(axis=0)
    
    df1, df2, df3 = dfs
    indexes_intersect = df1.index.intersection(df2.index).intersection(df3.index)
    print(indexes_intersect)
    corr_matrices = []
    for df in dfs:
        df = df.loc[indexes_intersect]
        # check here why it is empty
        corr_matrices.append(df.set_index('area').T.corr(method='pearson'))
    return corr_matrices
