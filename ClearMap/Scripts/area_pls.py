#!/usr/bin/env python3
import yaml


def main():
    import logging
    import coloredlogs
    import argparse
    import pandas as pd
    from utils_PLS import pls, bootstrap_test

    logger = logging.getLogger(__name__)
    logging.basicConfig(format='[%(funcName)s] - %(asctime)s - %(message)s', level=logging.INFO)
    coloredlogs.install(level='DEBUG', logger=logger)

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="input csv file", metavar='PATH')
    parser.add_argument('-o', '--output', help="output base path", metavar='PATH')
    parser.add_argument('-b', '--bootstrap', help="number of bootstrap iterations",
                        type=int, default=1000)
    parser.add_argument('-p', '--procrustes', help="enable Procrustes rotation in bootstrap",
                        action='store_true', default=False)
    parser.add_argument('-c', '--columns', help="number of descriptor columns",
                        type=int, default=2)
    args = parser.parse_args()

    logger.info('reading data...')
    data = pd.read_csv(args.input)  # csv file organized according to the example
    if 'Unnamed: 0' in data.columns:
        data.drop('Unnamed: 0', axis=1, inplace=True)  # remove index column if present

    # create dummy matrix containing the experimental group of each sample
    dum = pd.get_dummies(data, columns=['group'])
    n = data['group'].nunique()

    x = dum.iloc[:,args.columns:(data.shape[1]-1)]
    y = dum.iloc[:,-n:]
    xm = x.to_numpy()
    ym = y.to_numpy()

    logger.info('computing PLS...')
    u, s, v = pls(xm, ym)

    logger.info('PLS computed, now performing bootstrap...')
    vs = bootstrap_test(xm, ym, v, u, args.bootstrap, args.procrustes)

    regions = list(data.columns)
    regions = regions[3:]
    vpd = pd.DataFrame((v/vs), columns=regions)
    upd = pd.DataFrame(u, columns=y.columns)

    logger.info('saving output data...')
    vpd.to_csv(args.output + '_saliences.csv', index_label=False)
    upd.to_csv(args.output + '_contrasts.csv', index_label=False)


if __name__ == "__main__":
    main()
