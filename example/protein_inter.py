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
from tqdm import tqdm

#Silence MDTraj Warnings
warnings.filterwarnings("ignore")

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Protein Interactions')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-s', required=True, help= 'File containing sections of interest(txt)')
parser.add_argument('-n', required=True, help= 'Number of protein residues')

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
top = traj.topology
num_prot_res = args.n

#Load sections of interest file
sections = open(sect, 'r').readlines()

#Check if output should go to seperate directory
if not os.path.exists('./prot_inter/'):
    os.mkdir('prot_inter')

#Loop through sections of interest
for i in tqdm(range(len(sections))):
    #Input sections from file
    name1, name2, sect1, sect2 = load_data.read_sections(sections, i, miss_res, top, num_prot_res)
    res_pairs = list(product(sect1, sect2))
    if max(sect1) > traj.n_residues:
        print('Error Refernceing Atoms not in trajectory!')
        print(sect1)
        print(traj)
        exit()
    elif max(sect2) > traj.n_residues:
        print('Error Refernceing Atoms not in trajectory')
        print(traj)
        print(sect2)
        exit()

    #Compute smallest distance from heavy atoms
    [dist, pairs] = md.compute_contacts(traj, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)

    #Determine the % of the trajectory residues are in contact
    df_per = pd.DataFrame(np.nan, index = sect2, columns = sect1) 
    
    #Determine total interactions b/w sections per frame
    tot_inter = np.zeros(len(dist[:,0]))

    #Compute % contact
    ref_res = sect1[0]
    per_res, hc_res1, hc_res2, hc_per = [], [], [],[]
    for i in range(len(res_pairs)): #Loop through residue pairs
        dist_i = dist[:,i]
        contact = 0 #conter for protein--ligand contact
        for j in range(len(dist_i)): #Loop through time
            if dist_i[j] < 0.5:
                contact += 1
                tot_inter[j] += 1
        per_contact = 100 * (contact/len(dist_i))
        [res1, res2] = res_pairs[i]
        if res1 == ref_res:
            per_res.append(per_contact)
        else:
            df_per[ref_res] = per_res
            ref_res = res1
            per_res = []
            per_res.append(per_contact)
        if per_contact > 60 and abs(res2 - res1) > 2:
            hc_res1.append(res1)
            hc_res2.append(res2)
            hc_per.append(per_contact)
    
    #Reformat and save DF
    df_per.to_csv('prot_inter/' + name1 + '_' + name2 + '_inter_per.csv') 
    df_hc = pd.DataFrame({'Residue 1': hc_res1, 'Residue 2': hc_res2, 'Contact %': hc_per})
    df_hc.to_csv('prot_inter/' + name1 + '_' + name2 + '_high_contact.csv') 

    #Output total contacts for section
    df_tot = pd.DataFrame({'Total Interactions': tot_inter})
    df_tot.to_csv('prot_inter/' + name1 + '_' + name2 + '_tot_inter.csv') 

print('Protein Interaction Analysis Complete')
print('-------------------------------------------------------------')
