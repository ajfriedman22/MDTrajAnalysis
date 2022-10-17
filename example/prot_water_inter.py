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
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-l', required=False, type=int, default = 0, help= 'Ligand residue ID')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
lig = args.l

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_inter/')
import water_inter

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Test that ligand ID Assigned properly
test = traj.topology.select('resid ' + str(lig) + ' and resname LIG')
if test.size==0:
    print('Ligand not named correctly! Exiting Immediately!')
    exit()

#Determine water ligand contacts
if lig != 0:
    lig = lig
    water_indices, water_per_lig, water_neighbors = water_inter.water_neighbor(traj, lig, offset, False, 530)

    #Put residue indices in numerical order
    water_indices.sort()

    #Determine % contact with water residues
    water_res_per = np.zeros(len(water_indices))
    for j in range(len(water_neighbors)):#Loop through frames
        contact = water_neighbors[j]
        for n in range(len(contact)): #Loop through water interactions in the frame
            for i in range(len(water_indices)): #Loop through water residues to determine which is matched in frame
                if contact[n] == water_indices[i]:
                    water_res_per[i] += 1

    water_res_per = 100*(water_res_per/len(water_neighbors))

    #Output contacts
    output = open('Ligand_water_per.txt', 'w')
    output_high = open('Ligand_water_high.txt', 'w')
    for j in range(len(water_indices)):
        if water_res_per[j] > 10:
            output_high.write(str(water_indices[j]) + ': ' + str(water_res_per[j]) + '\n')
        output.write(str(water_indices[j]) + ': ' + str(water_res_per[j]) + '\n')
    output.close()
    output_high.close()
    print('Ligand Analysis Complete')

#Set residues of interest
res_interest = np.arange(1, lig-1, 1)
water_per_res = np.zeros(len(res_interest))
lig_res_per = np.zeros(len(res_interest))

#Loop through each residue of interest
for i in range(len(res_interest)):
    #Compute neighboring water molecules to residue of interest
    water_indices, water_per_res[i], lig_res_water_cont = water_inter.water_neighbor(traj, res_interest[i], offset, water_neighbors, 530)
    lig_res_per[i] = 100*sum(lig_res_water_cont)/len(lig_res_water_cont)

#Determine residues with highest number of water molecules interacting
max_water = max(water_per_res)
res_max_water = []
for i in range(len(res_interest)):
    if water_per_res[i] > 20:
        res_max_water.append(res_interest[i])
np.savetxt('res_max_water.txt', res_max_water)

#Ouput residues which simultaneously interact with the same water molecules as the ligand
output = open('lig_res_water.txt', 'w')
for i in range(len(lig_res_per)):
    if lig_res_per[i] > 25:
        output.write(str(res_interest[i]) + ': ' + str(lig_res_per[i]) + '\n')
np.savetxt('lig_res_water_per.txt', lig_res_per)
print('Protein Analysis Complete')
