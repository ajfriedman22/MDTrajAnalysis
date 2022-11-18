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
parser.add_argument('-sect', required=True, help= 'File containing sections of interest(txt)')
parser.add_argument('-f', required=False, type=int, default = 0, help= 'Frame for centroid output water molecule indices')

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
sect_name = args.sect
frame_cen = args.f

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_inter/')
import water_inter

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_ns = traj.remove_solvent()
print(traj)

if frame_cen != 0:
    cen_ind = []#Set aray for atom indices of water molecules for centroid

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine protein residues to examine based on sections of interest
input_sect = open(sect_name, 'r').readlines()
res_interest = []
for line in input_sect:
    line = line.split()
    sect_res = np.linspace(int(line[1]), int(line[2]), num = (int(line[2])-int(line[1])+1))
    #Add to final list
    for n in sect_res:
        res_interest.append(n)

#Put protein residue list array in numerical orde
res_interest.sort()

if lig != 0:
    #Test that ligand ID Assigned properly
    test = traj_ns.topology.select('resid ' + str(lig-offset) + ' and resname LIG')
    if test.size==0:
        print('Ligand not named correctly! Exiting Immediately!')
        exit()

    #Combine ligand and protein residues to make master list of residues of interest
    res_interest.append(lig)

#Loop through each residue of interest and determine water contacts
res_interest_water_neighbors = []
for i in range(len(res_interest)):
    #Compute neighboring water molecules to residue of interest
    water_neighbors = water_inter.water_neighbor(traj, res_interest[i], offset, lig)

    #Add to master list of water_neighbors
    res_interest_water_neighbors.append(water_neighbors)

#Determine number of frames in trajectory
frames = len(res_interest_water_neighbors[0])

#Determine percent water mediated contacts between all residues of interest
per_wat_contact = np.zeros((len(res_interest), len(res_interest)))
for i in range(len(res_interest)):
    water_contactA = res_interest_water_neighbors[i]
    for j in range(len(res_interest)):
        if i == j:
            per_wat_contact[i][j] = 100
        elif per_wat_contact[j][i] != 0:
            per_wat_contact[i][j] = per_wat_contact[j][i]
        else:
            water_contactB = res_interest_water_neighbors[j]
            #sort through frames to determine total contacts
            count = 0
            for t in range(frames):
                #determine if there are any common elements in this frame
                water_contactA_t = water_contactA[t]
                water_contactB_t = water_contactB[t]
                if np.in1d(water_contactA_t, water_contactB_t).any():
                    count += 1
                #Determine which water index is in contact in the frame of the centroid
                if t == frame_cen and frame_cen != 0 and j > i:
                    if np.in1d(water_contactA_t, water_contactB_t).any():
                        for n in water_contactA_t:
                            if n in water_contactB_t:
                                cen_ind.append(n)
                                break
                    else:
                        cen_ind.append(0)

            per_wat_contact[i][j] = np.round(100*count/frames, decimals = 2)
print(cen_ind)
print('Contacts Calculated')

#Print all present contacts to file
output = open('water_mediated_contact.txt', 'w')
output_cen = open('water_mediated_contact_cen_index.txt', 'w')
n = 0
for i in range(len(res_interest)):
    for j in range(len(res_interest)):
        if j > i:
            if per_wat_contact[i][j] > 25:
                output.write(str(res_interest[i]) + ' -- ' + str(res_interest[j]) + ': ' + str(per_wat_contact[i][j]) + '\n')
                output_cen.write(str(res_interest[i]) + ' -- ' + str(res_interest[j]) + ': ' + str(cen_ind[n]) + '\n')
            n += 1
print('Output files written')
