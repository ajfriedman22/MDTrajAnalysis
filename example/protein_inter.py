#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, protein and ligand RMSD, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-sect', required=True, help= 'File containing sections of interest(txt)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
sect = args.sect
if sect.split('.')[-1] != 'txt': #Add default file extension if not in input
    sect = sect + '.txt'

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_ns = traj.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
top = traj_ns.topology
del traj

#Load sections of interest file
sections = open(sect, 'r').readlines()

#Loop through sections of interest
for i in range(len(sections)):
    [name1, sect1_start, sect1_end, name2, sect2_start, sect2_end] = sections[i].split()
    
    #Compute distance between all residues in sect1 and sect2
    sect1_start = int(sect1_start)-1-miss_res
    sect1_end = int(sect1_end)-1-miss_res
    sect2_start = int(sect2_start)-1-miss_res
    sect2_end = int(sect2_end)-1-miss_res

    sect1 = np.linspace(sect1_start, sect1_end, num=sect1_end-sect1_start+1)
    sect2 = np.linspace(sect2_start, sect2_end, num=sect2_end-sect2_start+1)
    res_pairs = list(product(sect1, sect2))
    [dist, pairs] = md.compute_contacts(traj_ns, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)
    
    #Determine the % of the trajectory residues are in contact
    output_per = open(name1 + '_' + name2 + '_inter_per.txt', 'w')
    output_high_contact = open(name1 + '_' + name2 + '_high_contact.txt', 'w')
    
    #Determine total interactions b/w sections per frame
    tot_inter = np.zeros(len(dist[:,0]))

    #Compute % contact
    for i in range(len(res_pairs)):
        dist_i = dist[:,i]
    
        contact = 0 #conter for protein--ligand contact
        for j in range(len(dist_i)):
            if dist_i[j] < 0.4:
                contact += 1
                tot_inter[j] += 1
        per_contact = 100 * (contact/len(dist_i))
        [res1, res2] = res_pairs[i]
        output_per.write(str(res1-1-miss_res) + ' -- ' + str(res2-1-miss_res) + ': ' + str(per_contact)+ '\n')
        if per_contact > 75:
            output_high_contact.write(str(res_pairs[i]) + '\n')
    
    #Output total contacts for section
    np.savetxt(name1 + '_' + name2 + '_tot_inter.txt', tot_inter)
