#!/ usr / bin / env python
import numpy as np
import argparse
import sys
import os.path
import time
import pandas as pd
import mdtraj as md
from itertools import product
from tqdm import tqdm

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
traj = load_data.mdtraj_load(File_traj, File_gro, False, False)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine protein residues to examine based on sections of interest
input_sect = open(sect_name, 'r').readlines()
res_interest, compute_lig = [],[]
for i, line in enumerate(input_sect):
    sect_res = line.split()
    #Add to final list
    for n in sect_res:
        res_interest.append(int(n))
        if i == 0:
            compute_lig.append(True)
        else:
            compute_lig.append(False)

if lig != 0:
    #Test that ligand ID Assigned properly
    lig_atoms = traj.topology.select('resid ' + str(lig-offset) + ' and resname ' + lig_name)
    if lig_atoms.size==0:
        raise Exception('Ligand not named correctly! Exiting Immediately!')
    
    water_neighbors_lig = water_inter.water_neighbor(traj, lig, offset, lig)
load_time = time.time()

#Loop through each residue of interest and determine water contacts
res_interest_water_neighbors = []
for i in range(len(res_interest)):
    #Compute neighboring water molecules to residue of interest
    water_neighbors = water_inter.water_neighbor(traj, res_interest[i], offset, lig)

    #Add to master list of water_neighbors
    res_interest_water_neighbors.append(water_neighbors)

if len(res_interest_water_neighbors) == 0:
    raise Exception('Error No Water Neighbors Found! Exiting Immediately!')
neighbor_time = time.time()
print('Water Neighbors Found')

main_df = pd.DataFrame()
for r, res in enumerate(res_interest):
    water_neighbors_res = res_interest_water_neighbors[r]
    res_atoms = traj.topology.select(f'resid {res-offset}')
    if compute_lig[r]:
        dist_tot_lig = np.zeros(traj.n_frames)
    dist_tot_res = np.zeros(traj.n_frames)
    for t in tqdm(range(traj.n_frames)):
        #Get just this frame from the trajectory
        traj_t = traj.slice(t)
        if compute_lig[r]:
            #Determine common water neighbors during this frame
            water_neighbors_lig_t = water_neighbors_lig[t]
            water_neighbors_res_t = water_neighbors_res[t]
            common_neighbors = set(water_neighbors_lig_t) & set(water_neighbors_res_t)
            if len(common_neighbors) != 0:
                #Get indices for water molecules
                water_index = []
                for wat in common_neighbors:
                    wat_atoms = traj.topology.select(f'resid {wat-offset}')
                    for atom in wat_atoms:
                        water_index.append(atom)
                #Determine minimum distance between water molecule and ligand
                lig_pairs = list(product(lig_atoms, water_index))
                dist_lig = md.compute_distances(traj_t, atom_pairs=lig_pairs)
                dist_tot_lig[t] = np.min(dist_lig[0])
            
                #Determine minimum distance between water molecule and protein residue
                res_pairs = list(product(res_atoms, water_index))
                dist_res = md.compute_distances(traj_t, atom_pairs=res_pairs)
                dist_tot_res[t] = np.min(dist_res[0])
            else:
                dist_tot_lig[t] = 0.4
                dist_tot_res[t] = 0.4
        else:
            #Get water atoms
            water_neighbors_res_t = water_neighbors_res[t]
            for wat in water_neighbors_res_t:
                wat_atoms = traj.topology.select(f'resid {wat-offset}')
                for atom in wat_atoms:
                    water_index.append(atom)
            #Determine minimum distance between water molecule and protein residue
            res_pairs = list(product(res_atoms, water_index))
            dist_res = md.compute_distances(traj_t, atom_pairs=res_pairs)
            dist_tot_res[t] = np.min(dist_res[0])

    if compute_lig[r]:
        df = pd.DataFrame({'Distance': dist_tot_lig*10, 'Side': ['Ligand']*len(dist_tot_lig), 'Residue': [res]*len(dist_tot_lig)})
        main_df = pd.concat([main_df, df])
    df2 = pd.DataFrame({'Distance': dist_tot_res*10, 'Side': ['Residue']*len(dist_tot_res), 'Residue': [res]*len(dist_tot_res)})
    main_df = pd.concat([main_df, df2])
main_df.to_csv('water_med_inter_dist.csv')
print(main_df.groupby('Residue')['Distance'].mean())