#!/usr/bin/env python3
import yaml

##################################################################
# Core Area PLS developed originally by Ludovico Silvestri
# Documentation by AS
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
