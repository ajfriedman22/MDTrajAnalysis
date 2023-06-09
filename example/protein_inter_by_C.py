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
parser.add_argument('-f', required=True, help= 'input file (line 1: primary res name_1 name_2 ... line 2: res_interest_1 res_interest_2 ...)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
miss_res = args.m
sect = args.f
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
input_sect = open(sect, 'r').readlines()
main_res = int(input_sect[0].split(' ')[0]) - miss_res - 1
main_name = input_sect[0].split(' ')[1:]
res_interest = input_sect[1].split(' ')

#Determine atom numbers for main residue
main_atom = []
for i in main_name:
    atom = top.select('resid ' + str(main_res) + ' and name ' + i)
    if len(atom) == 0:
        print('Residue ' + str(main_res) + ' has no atom named ' + str(i))
    else:
        main_atom.append(atom[0])

#Declare empty vector for percent interactions are formed
inter_by_C = np.zeros((len(main_atom), len(res_interest)))

#Check if output should go to seperate directory
if os.path.exists('./prot_inter/'):
    dir_name = 'prot_inter/'
else:
    dir_name = ''

progress = ProgressBar(len(res_interest), fmt=ProgressBar.FULL)

#Loop through sections of interest
for i in range(progress.total):
    #Determine atom numbers for residue or interest
    interest_atom = top.select("resid " + str(int(res_interest[i])-1-miss_res) + " and mass > 4.0")
    if len(interest_atom) == 0:
        print('No heacy atoms in residue ' + str(res_interest[i]))
        break
    
    #Create atom pairs
    res_pairs = list(product(main_atom, interest_atom))

    #Compute smallest distance from heavy atoms
    dist = md.compute_distances(traj, res_pairs)
    time, num_pairs = np.shape(dist)
    
    #Compute % contact
    per_res = []
    for m in range(len(main_atom)): #Loop through pairs
        count = 0
        for t in range(time):
            for j in range(len(interest_atom)):
                n = int(m*len(interest_atom) + j)
                if dist[t][n] < 0.5:
                    count += 1
                    break
        per_res.append(100*count/time)
    inter_by_C[:,i] = per_res
    progress.current += 1
    progress()

#Save to df
df = pd.DataFrame(inter_by_C, index = main_name, columns = res_interest)
df.to_csv(dir_name + 'inter_by_C.csv')
progress.done()
print('Protein Interaction Analysis Complete')
print('-------------------------------------------------------------')
