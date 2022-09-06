#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from sklearn.decomposition import PCA
from itertools import combinations
import argparse
from itertools import product
from statistics import stdev
import sys
from mpl_toolkits import mplot3d
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, protein and ligand RMSD, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, default = 0, help= 'Supply the number of missing terminal residues')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m

#Source custom functions
sys.path.insert(1, '$PROJECT/code/MDTrajAnalysis/Traj_process/')
import traj 

sys.path.insert(1, '$PROJECT/code/MDTrajAnalysis/protein_inter/')
import water_inter

#Load Trajectory
traj = traj.mdtraj_load(File_traj, File_gro)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine the distance b/w all water and protein residues
water_dist = water_inter.prot_water_dist(traj, [0,180], offset)

print(np.shape(water_dist))
