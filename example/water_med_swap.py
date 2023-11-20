#!/ usr / bin / env python
import numpy as np
import argparse
import sys
import os.path
import time
import pandas as pd
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
parser.add_argument('-sf', required=False, default=False, help='Is input in section format? (Otherwise list individual residues)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
miss_res = args.m
lig = args.l
lig_name = args.ln
sect_name = args.s
sect_format = args.sf

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
    if sect_format:
        sect_res = np.linspace(int(line[0]), int(line[1]), num = (int(line[1])-int(line[0])+1))
        #Add to final list
        for n in sect_res:
            res_interest.append(int(n))
    else:
        for n in line:
            res_interest.append(int(n))

#Put protein residue list array in numerical order
res_interest.sort()

#Test that ligand ID Assigned properly
test = traj.topology.select('resid ' + str(lig-offset) + ' and resname ' + lig_name)
if test.size==0:
    raise Exception('Ligand not named correctly! Exiting Immediately!')

#Get water mediated interactions with ligand
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

#Determine number of frames in trajectory
frames = len(res_interest_water_neighbors[0])

#Determine percent water mediated contacts between all residues of interest
transition_time = []
for i in tqdm(range(len(res_interest))):
    trans_time_res = []
    water_contactA = res_interest_water_neighbors[i]
            
    #sort through frames to determine frames per swap
    common_contact = []
    #Get common waters
    for t in range(frames):
        #determine common elements in this frame
        water_contactA_t = water_contactA[t]
        water_contact_lig_t = water_neighbors_lig[t]
        prot_water_lig = set(water_contactA_t) & set(water_contact_lig_t)
        common_contact.append(prot_water_lig)
    #set buffer for frames
    buffer = 5
    count = 0
    ref_water = common_contact[0]
    for t in range(frames-buffer):
        for b in range(buffer):
            if set(ref_water) & set(common_contact[t+b]):
                count += 1+b
                t=t+b
                break
            elif b==buffer:
                if count > 0:
                    trans_time_res.append(count/frames*300)
                count = 0
                ref_water = common_contact[t]
    while 0.0 in trans_time_res:
        trans_time_res.remove(0.0)
    transition_time.append(trans_time_res)
contact_time = time.time()
print('Contacts Calculated')

#Print all present contacts to file

main_df = pd.DataFrame()
for r, res in enumerate(res_interest):
    df = pd.DataFrame({'Residue': [res]*len(transition_time[r]), 'Transition time': transition_time[r]})
    main_df = pd.concat([main_df, df])
main_df.to_csv('water_med_trans_time.csv')
print(main_df.groupby('Residue')['Transition time'].mean())
output_time = time.time()
print('Output files written')

print('Load data: ' + str(load_time - start_time))
print('Compute Contacts: ' + str(neighbor_time - load_time))
print('Determine Contacts: ' + str(contact_time - neighbor_time))
print('Output data: ' + str(output_time - contact_time))
print('-----------------------------------------------------------------')

