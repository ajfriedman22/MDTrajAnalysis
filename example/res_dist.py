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
from tqdm import tqdm

#Silence MDTraj Warnings
warnings.filterwarnings("ignore")

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Protein Interactions')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-s', required=True, help= 'File containing sections of interest(txt)')

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
offset = 1 + miss_res

#Load sections of interest file
sections = open(sect, 'r').readlines()
res_pairs = np.zeros((len(sections), 2))
for i, line in enumerate(sections):
    resid1, atom1, resid2, atom2 = line.split(' ')
    res1 = traj.topology.select(f'resid {int(resid1)-offset} and name {atom1.strip()}')
    res2 = traj.topology.select(f'resid {int(resid2)-offset} and name {atom2.strip()}')
    if len(res1) != 1 and len(res2) != 1:
        raise Exception(f'')
    res_pairs[i,:] = [int(res1[0]), int(res2[0])]

#Determine the distance 
dist = md.compute_distances(traj, atom_pairs=res_pairs)
dist = dist*10

#Save to file
main_df = pd.DataFrame()
for i, line in enumerate(sections):
    resid1, atom1, resid2, atom2 = line.split(' ')
    dist_i = dist[:,i]
    print(f'{resid1}-{resid2}: {np.mean(dist_i)}')
    df = pd.DataFrame({'Distance(A)': dist_i, 'Resid 1': [resid1]*len(dist_i), 'Resid 2': [resid2]*len(dist_i)})
    main_df = pd.concat([main_df, df])
main_df.to_csv(f'dist_lig_coord_water.csv')