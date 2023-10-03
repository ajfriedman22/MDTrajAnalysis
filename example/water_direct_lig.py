#!/ usr / bin / env python
import numpy as np
import argparse
import sys
import os.path
import time
from itertools import product
import pandas as pd
import mdtraj as md

start_time = time.time()

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Water Mediated Interactions')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-l', required=False, type=int, default = 0, help= 'Ligand residue ID')
parser.add_argument('-ln', required=False, type=str, default = 'LIG', help= 'Ligand name')
parser.add_argument('-s', required=True, help= 'File containing sections of interest(txt)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
miss_res = args.m
lig = args.l
lig_name = args.ln
sect_name = args.s

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_analysis/')
import water_inter

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro, False, True)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine protein residues to examine based on sections of interest
input_sect = open(sect_name, 'r').readlines()
res_interest = []
for line in input_sect:
    line = line.split()
    if len(line) > 2:
        res_interest_wat = line
    else:
        sect_res = np.linspace(int(line[0]), int(line[1]), num = (int(line[1])-int(line[0])+1))
        #Add to final list
        for n in sect_res:
            res_interest.append(int(n))

#Put protein residue list array in numerical order
res_interest.sort()

#Get all heavy atoms for residues of interest
prot_heavy = []
for res in res_interest:
    atoms = traj.topology.select(f'resid {res-offset} and element != H')
    for atom in atoms:
        prot_heavy.append(atom)

#Test that ligand ID Assigned properly
test = traj.topology.select('resid ' + str(lig-offset) + ' and resname ' + lig_name)
if test.size==0:
    raise Exception('Ligand not named correctly! Exiting Immediately!')
#Get all atoms engaged in direct non-bonded interacts with the ligand
lig_atom = traj.topology.select(f'resid {lig-offset} and element != H')

#Get water mediated interactions with ligand
water_neighbors_lig = (water_inter.water_neighbor(traj, lig, offset, lig))

#Loop through each residue of interest and determine water contacts
res_interest_water_neighbors = []
for i in range(len(res_interest_wat)):
    #Compute neighboring water molecules to residue of interest
    water_neighbors = water_inter.water_neighbor(traj, float(res_interest_wat[i]), offset, lig)

    #Add to master list of water_neighbors
    res_interest_water_neighbors.append(water_neighbors)

if len(res_interest_water_neighbors) == 0:
    raise Exception('Error No Water Neighbors Found! Exiting Immediately!')
print('Water Neighbors Found')

#Determine number of frames in trajectory
frames = len(res_interest_water_neighbors[0])

#Compute distances to all protein residues
atom_pairs = list(product(lig_atom, prot_heavy))
dist = md.compute_distances(traj, atom_pairs, periodic=True)

#Get all atoms 
direct_contact_atom = [] #list for atoms in direct contact with ligand
for t in range(frames):
    atom_contact = []
    for i in range(len(atom_pairs)):
        if dist[t,i] < 0.5 and atom_pairs[i][1] not in atom_contact:
            atom_contact.append(atom_pairs[i][1])
    direct_contact_atom.append(atom_contact)
print('Direct Interacting Atoms Found')

#Determine percent of direct and water mediated interactions formed with each ligand heavy atom
for t in range(frames):
    traj_frame = traj.slice(t)
    for i in range(len(res_interest_wat)):
        if t==0:
            globals()[f'atom_water_{res_interest_wat[i]}'] = []
        #Seperate water contacts for this residue
        water_contact_res = res_interest_water_neighbors[i]
        water_common = set(water_neighbors_lig[t]) & set(water_contact_res[t])
        if set(water_neighbors_lig[t]) & set(water_contact_res[t]):
            #Get minimum heavy atom distance for residue of interest
            res_heavy = traj.topology.select(f'resid {res_interest[i]-offset} and element != H')
            water_atom = []
            for wat in water_common:
                water_atom.append(traj.topology.select(f'resid {wat-offset}')[0])
            atom_pairs = list(product(res_heavy, water_atom))
            dist = md.compute_distances(traj_frame, atom_pairs, periodic=True)
            globals()[f'atom_water_{res_interest_wat[i]}'].append(atom_pairs[np.argmin(dist[0,:])][0])
print('Water Mediated Atoms Found')

#Determine minimum distance between residues engaged in water mediated and direct intereactions
df = pd.DataFrame()
for res in res_interest_wat:
    globals()[f'min_neighbor_dist_{res}'] = []
    for t in range(len(globals()[f'atom_water_{res}'])):
        atom_pairs = list(product(direct_contact_atom[t], [globals()[f'atom_water_{res}'][t]]))
        traj_frame = traj.slice(t)
        dist = md.compute_distances(traj_frame, atom_pairs)
        globals()[f'min_neighbor_dist_{res}'].append(np.min(dist[0]))
    df[res] = globals()[f'min_neighbor_dist_{res}']

#Print all present contacts to file
df.to_csv('water_direct_lig.csv')

print('Output files written')

