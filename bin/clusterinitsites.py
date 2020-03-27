#!/usr/bin/env python

import pandas as pd
import argparse as ap
import subprocess as sp
import ntpath
import logging
import os

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = ap.ArgumentParser()
parser.add_argument('-p', '--peakfile', required = True,
                    help = 'narrowPeak file containing initiation site coordinates')
parser.add_argument('-f', '--format', default = 'narrow', choices = ['narrow', 'BED4', 'BED6'],
                    help = 'format of the given peakfile')
parser.add_argument('-n', '--sitespercluster', default = 2, type = int,
                    help = 'minimum number of sites per identified cluster')
parser.add_argument('-o', '--outdir', help = 'directory to which results should be written', default = '.')
args = parser.parse_args()

colnames = {'narrow': ['chr', 'start', 'end', 'name', 'score', 'strand', 'fold', 'logp', 'logq', 'relsummittostartpos'],
            'BED4': ['chr', 'start', 'end', 'name'],
            'BED6': ['chr', 'start', 'end', 'name', 'score', 'strand']}

suffix = '.narrowPeak' if args.format == 'narrow' else '.bed'

# reading in MACS output
logging.info('generating required annotation file')
bed = pd.read_csv(args.peakfile, sep = '\t', header = None)
bed.columns = colnames[args.format]
bed['featurename'] = 'ispeak'

# generating required annotation file
with open(args.peakfile.replace(suffix, '.csa'), 'w') as csa:
    csa.write('\t'.join(['peakid', 'featuretype']) + '\n')
    bed[['name', 'featurename']].to_csv(csa, sep = '\t', index = False, mode = 'a', header = False)

# calculating pairwise distances between MACs peaks
logging.info('computing median distance')
distances = []
for group in bed.groupby('chr'):
    frame = group[1]
    for i in frame.index[:-1]:
        distances.append(frame.at[i + 1, 'start'] - frame.at[i, 'end'])

median = pd.Series(distances).median()
logging.info('median distance between peaks is %d' % median)

# clustering initiation sites by distance
logging.info('calling clusterscan.py clusterdist %s %s -a %s -d %i -n %i' %
             (args.peakfile, args.peakfile.replace(suffix, '.csa'),
              os.path.join(args.outdir, '_'.join(ntpath.basename(args.peakfile).split('_')[:-1])),
              median, args.sitespercluster))
subprocess = sp.Popen('clusterscan.py clusterdist %s %s -a %s -d %i -n %i' %
                      (args.peakfile, args.peakfile.replace(suffix, '.csa'),
                       os.path.join(args.outdir, '_'.join(ntpath.basename(args.peakfile).split('_')[:-1])),
                       median, args.sitespercluster),
                      shell = True)
subprocess.wait()

# processing results
resultbed = pd.read_csv(os.path.join(args.outdir,
                                     '_'.join(ntpath.basename(args.peakfile).split('_')[:-1] + ['clusters.bed'])),
                        sep = '\t', header = None, skiprows = 1)
resultbed.columns = ['chr', 'start', 'end', 'name', 'nsites', 'strand', 'featurename']
resultbed.loc[:, 'name'] = ['_'.join([args.peakfile.split('_')[2], str(i)]) for i in range(1, len(resultbed) + 1)]
resultbed.to_csv(os.path.join(args.outdir, '_'.join(ntpath.basename(args.peakfile).split('_')[:-1] + ['clusters.bed'])),
                 header = None, index = None, sep = '\t')
