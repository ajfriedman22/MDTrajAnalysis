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
parser.add_argument('-m', required=False, default = 0, help= 'Supply the number of missing terminal residues(default 0)')

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
prefix = '/ocean/projects/cts160011p/afriedma/code/MDTrajAnalysis/'
sys.path.insert(1, prefix + 'Traj_process/')
import traj 

sys.path.insert(1, prefix + 'protein_inter/')
import water_inter

#Load Trajectory
traj = traj.mdtraj_load(File_traj, File_gro)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Set residues of interest
res_interest = [121]

#Loop through each residue of interest
for i in res_interest:
    #Compute neighboring water molecules to residue of interest
    water_indices = water_inter.water_neighbor(traj, 121, offset)

    #Determine the distance b/w all water and protein residues
    water_dist = water_inter.prot_water_dist(traj, 121, offset, water_indices)

    #Determine Number of H-bonds and VDW interactions b/w residue and water molecules in each frame


    print(np.shape(water_dist))
