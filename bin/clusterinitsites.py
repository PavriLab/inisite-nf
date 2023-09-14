#!/usr/bin/env python

import pandas as pd
import argparse as ap
import subprocess as sp
import logging
import os

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = ap.ArgumentParser()
parser.add_argument('-p', '--peakfile', required = True,
                    help = 'narrowPeak file containing initiation site coordinates')
parser.add_argument('-n', '--sitespercluster', default = 2, type = int,
                    help = 'minimum number of sites per identified cluster')
parser.add_argument('--clusterScan', required = True,
                    help = 'path to clusterscan.py script')
parser.add_argument('-o', '--outputPrefix', required = True,
                    help = 'prefix for the output file')
args = parser.parse_args()

# reading in MACS output
logging.info('generating required annotation file')
bed = pd.read_csv(args.peakfile, sep = '\t', header = None, usecols = [0, 1, 2, 3],
                  names = ['chr', 'start', 'end', 'name'])
bed['featurename'] = 'ispeak'

# generating required annotation file
annotationFile = args.peakfile.split('.')[0] + '.csa'
with open(annotationFile, 'w') as csa:
    csa.write('\t'.join(['peakid', 'featuretype']) + '\n')
    bed[['name', 'featurename']].to_csv(csa, sep = '\t', index = False, mode = 'a', header = False)

# calculating pairwise distances between MACS peaks
logging.info('computing median distance')
distances = []
for _, group in bed.groupby('chr'):
    for i in group.index[:-1]:
        distances.append(group.at[i + 1, 'start'] - group.at[i, 'end'])

median = int(pd.Series(distances).median())
logging.info('median distance between peaks is %d' % median)

logging.info('clustering with clusterscan.py clusterdist {0} {1} -a {2} -d {3} -n {4}'.format(
                            args.peakfile,
                            annotationFile,
                            args.outputPrefix,
                            median,
                            args.sitespercluster))

subprocess = sp.Popen('{0} clusterdist {1} {2} -a {3} -d {4} -n {5}'.format(
                            args.clusterScan,
                            args.peakfile,
                            annotationFile,
                            args.outputPrefix,
                            median,
                            args.sitespercluster),
                      shell = True)
subprocess.wait()

# processing results
resultbed = pd.read_csv(args.outputPrefix + '_clusters.bed', sep = '\t', header = None, skiprows = 1)
resultbed.columns = ['chr', 'start', 'end', 'name', 'nsites', 'strand', 'featurename']
resultbed.loc[:, 'name'] = ['iniZone_{0}'.format(i) for i in range(1, len(resultbed) + 1)]
resultbed.to_csv(args.outputPrefix + '_clusters.bed', header = None, index = None, sep = '\t')
