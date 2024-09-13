#!/usr/bin/env python3
import yaml

##################################################################
# PLS Analysis and Bootstrap Testing Script

# This Python script performs Partial Least Squares (PLS) analysis and bootstrap testing on experimental data provided in CSV format. It processes input data, conducts statistical analysis, and saves the results to specified output files.

# ## Overview

# The script is designed to:
# - Read data from a CSV file.
# - Perform Partial Least Squares (PLS) analysis to identify relationships between input variables (X) and response variables (Y).
# - Conduct bootstrap testing to assess the stability and reliability of the PLS results.
# - Optionally apply Procrustes rotation during bootstrap iterations to refine the analysis.
# - Output the saliences and contrasts resulting from the PLS analysis as CSV files.

# ## Key Features

# - **Command-line Interface**: Uses command-line arguments to specify input files, output paths, and analysis parameters.
# - **Data Preprocessing**: Automatically handles index columns and creates dummy matrices for categorical data such as experimental groups.
# - **PLS Analysis**: Decomposes the relationship between input data (X) and output data (Y) using PLS to generate matrices describing these relationships.
# - **Bootstrap Testing**: Repeats the analysis over multiple iterations to test the stability of the results. Procrustes rotation can be applied during this process.
# - **Results Export**: Outputs saliences and contrasts to CSV files for further analysis.

# ## Parameters

# The script accepts the following command-line arguments:

# - `-i`, `--input`: Path to the input CSV file containing the data.
# - `-o`, `--output`: Base path for the output files. Results will be saved with this prefix.
# - `-b`, `--bootstrap`: Number of bootstrap iterations (default is 5000).
# - `-p`, `--procrustes`: Enables Procrustes rotation during bootstrap testing (optional, default is False).
# - `-c`, `--columns`: Number of descriptor columns to use in the PLS analysis (default is 2).

# ## Workflow

# 1. **Input Data**: Reads a CSV file containing experimental data, assuming one column contains group labels and others contain descriptors or measurements.
# 2. **Data Processing**: Converts categorical group data into a dummy matrix using one-hot encoding and selects appropriate columns for PLS analysis.
# 3. **PLS Analysis**: Applies the PLS method to identify relationships between the independent variables (X) and dependent variables (Y), generating matrices `U`, `S`, and `V`.
# 4. **Bootstrap Testing**: Conducts a bootstrap analysis with multiple iterations (default 5000) to test the stability of the `V` matrix from PLS. Optionally applies Procrustes rotation.
# 5. **Result Export**: Outputs salience values and contrasts to CSV files for further inspection.

# ## Example Usage

# ```bash
# python3 script.py -i data/input.csv -o results/output -b 10000 -p -c 3

##################################################################

def main():
    
    import logging
    import coloredlogs
    import argparse
    import pandas as pd
    from utils_PLS import pls, bootstrap_test

    """
    Main function to perform Partial Least Squares (PLS) analysis and bootstrap testing on input data.

    This function reads input data from a CSV file, processes it, performs PLS analysis, and then
    conducts bootstrap testing to assess the stability of the PLS results. The results are saved to
    specified output files.

    Parameters:
    None (the function uses command-line arguments for input parameters)

    Returns:
    None
    """
    
    #logger = logging.getLogger(__name__)
    #logging.basicConfig(format='[%(funcName)s] - %(asctime)s - %(message)s', level=logging.INFO)
    #coloredlogs.install(level='DEBUG', logger=logger)

    # Setup argument parser for command-line options, arguments, and sub-commands
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', 
                        help="input csv file", metavar='PATH')
    parser.add_argument('-o', '--output', 
                        help="output base path", metavar='PATH')
    parser.add_argument('-b', '--bootstrap', 
                        help="number of bootstrap iterations", type=int, default=5000)
    parser.add_argument('-p', '--procrustes', 
                        help="enable Procrustes rotation in bootstrap", action='store_true', default=False)
    parser.add_argument('-c', '--columns', 
                        help="number of descriptor columns", type=int, default=2)
    args = parser.parse_args()

    # Read input data from CSV file
    data = pd.read_csv(args.input)
    
    # Drop the 'Unnamed: 0' column if it exists (often an index column from pandas)
    if 'Unnamed: 0' in data.columns:
        data.drop('Unnamed: 0', axis=1, inplace=True)

    # Create a dummy matrix for the experimental groups
    dum = pd.get_dummies(data, columns=['group'])
    n = data['group'].nunique()  # Number of unique groups

    # Select columns for the PLS analysis
    x = dum.iloc[:, args.columns:(data.shape[1] - 1)]
    y = dum.iloc[:, -n:]
    xm = x.to_numpy()
    ym = y.to_numpy()

    # Perform PLS analysis
    u, s, v = pls(xm, ym)

    # Perform bootstrap testing
    vs = bootstrap_test(xm, ym, v, u, args.bootstrap, args.procrustes)

    # Define regions based on data columns
    regions = list(data.columns)
    regions = regions[3:]  # Adjust to match the specific columns of interest

    # Create DataFrames for PLS results and bootstrap results
    vpd = pd.DataFrame((v / vs), columns=regions)
    upd = pd.DataFrame(u, columns=y.columns)

    # Save the results to CSV files
    vpd.to_csv(args.output + '_saliences.csv', index=False)
    upd.to_csv(args.output + '_contrasts.csv', index=False)


if __name__ == "__main__":
    
    import random
    random.seed(0)
    import numpy as np
    np.random.seed(0)
    main()
