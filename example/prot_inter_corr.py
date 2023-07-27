#!/ usr / bin / env python
from __future__ import print_function
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product
import warnings
import pandas as pd
import re
from scipy.stats import chi2_contingency

class ProgressBar(object):
    DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
    FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'

    def __init__(self, total, width=40, fmt=DEFAULT, symbol='=',
                 output=sys.stderr):
        assert len(symbol) == 1

        self.total = total
        self.width = width
        self.symbol = symbol
        self.output = output
        self.fmt = re.sub(r'(?P<name>%\(.+?\))d',
            r'\g<name>%dd' % len(str(total)), fmt)

        self.current = 0

    def __call__(self):
        percent = self.current / float(self.total)
        size = int(self.width * percent)
        remaining = self.total - self.current
        bar = '[' + self.symbol * size + ' ' * (self.width - size) + ']'

        args = {
            'total': self.total,
            'bar': bar,
            'current': self.current,
            'percent': percent * 100,
            'remaining': remaining
        }
        print('\r' + self.fmt % args, file=self.output, end='')

    def done(self):
        self.current = self.total
        self()
        print('', file=self.output)

#Silence MDTraj Warnings
warnings.filterwarnings("ignore")

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, protein and ligand RMSD, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-s', required=True, help= 'File containing residues to determine correlated presence of interactions(txt)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
miss_res = args.m
sect = args.s
if sect.split('.')[-1] != 'txt': #Add default file extension if not in input
    sect = sect + '.txt'

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro, True, True)
top = traj.topology

#Load sections of interest file
sections = open(sect, 'r').readlines()

#Check if output should go to seperate directory
if os.path.exists('./prot_inter/'):
    dir_name = 'prot_inter/'
else:
    dir_name = ''

progress = ProgressBar(len(sections), fmt=ProgressBar.FULL)

#Loop through sections of interest
for i in range(progress.total):
    #Input sections from file
    if len(sections[i].split()) != 3:
        raise Exception('Invalid Input File Format Expected three inputs per line not ' + len(sections[i].split()))
    else:
        main_inter, target1, target2 = sections[i].split()
    res_pairs = [[int(main_inter)-1-miss_res, int(target1)-1-miss_res],[int(main_inter)-1-miss_res, int(target2)-1-miss_res]]

    #Compute smallest distance from heavy atoms
    [dist, pairs] = md.compute_contacts(traj, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)
    time, pairs = np.shape(dist)

    #Compute % contact
    inter = np.zeros((2,2))
    for t in range(time):
        if dist[t][0] <= 0.5 and dist[t][1] <= 0.5:
            inter[1][1] += 1
        elif dist[t][0] <= 0.5 and dist[t][1] > 0.5:
            inter[1][0] += 1
        elif dist[t][0] > 0.5 and dist[t][1] <= 0.5:
            inter[0][1] += 1
        else:
            inter[0][0] += 1

    print(inter)
    #Determine interaction correlation coefficent
    chi = chi2_contingency(inter)
    print(chi[1])

    #Update Progress
    progress.current += 1
    progress()
progress.done()
print('Protein Interaction Analysis Complete')
print('-------------------------------------------------------------')
