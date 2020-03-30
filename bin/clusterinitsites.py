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
parser.add_argument('-o', '--outputPrefix', required = True,
                    help = 'prefix for the output file')
args = parser.parse_args()

# reading in MACS output
logging.info('generating required annotation file')
bed = pd.read_csv(args.peakfile, sep = '\t', header = None, usecols = [0, 1, 2, 4],
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
for group in bed.groupby('chr'):
    frame = group[1]
    for i in frame.index[:-1]:
        distances.append(frame.at[i + 1, 'start'] - frame.at[i, 'end'])

median = pd.Series(distances).median()
logging.info('median distance between peaks is %d' % median)

logging.info('clustering initiation sites')
subprocess = sp.Popen('clusterscan.py clusterdist {0} {1} -a {2} -d {3} -n {4}'.format(
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
