#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
import pandas as pd
import warnings

#Silence MDTraj Warnings
warnings.filterwarnings("ignore")

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of SASA for different protein sections')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-s', required=False, default = 'none', type = str, help= 'File name for residue ranges for computed total SASA (name initial final)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
sect = args.s

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
top = traj.topology

#Set protein offset based on missing residues
offset = 1 + miss_res

#Seperate section for DSSP calculation
input_file = open(sect, 'r').readlines()

#Compute SASA for sections
df = pd.DataFrame()
for i in range(len(input_file)):
    if len(input_file[i].split()) == 3:
        [name, first, last] = input_file[i].split()
        exclude = None
    else:
        line = input_file[i].split()
        [name, first, last] = input_file[i].split()[0:3]
        exclude = np.array(line[4:], dtype=int)
        
    first, last = int(first), int(last)
    
    #Seperate indices for atoms of interest
    ind = np.linspace(first, last, num=last-first+1, dtype=int)
    if exclude != None:
        ind = np.setdiff1d(ind,exclude)
    traj_sect = traj.atom_slice(ind)
    
    #Compute SASA for section
    sasa_sect = md.shrake_rupley(traj_sect) 

    total_sasa = sasa_sect.sum(axis=1)
    df[name] = total_sasa
df.to_csv('sasa.csv')
