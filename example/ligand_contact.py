#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, protein and ligand RMSD, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-l', required=False, default = 'none', type = str, help= 'Ligand name if ligand analysis should be performed')

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

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
top = traj.topology


#Determine indices of protein and ligand residues in topology
prot_res = traj_ns.topology.select('protein')
lig_res = traj_ns.topology.select('resname ' + lig) #Select the ligand by name (based on topology) from the trajectory

#Determine distance between all protein residues and ligand
res_pairs = list(product(lig_res, prot_res))
[dist, pairs] = md.compute_contacts(traj_ns, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)

#Output % contact for each residue
output_per = open('residue_lig_per.txt', 'w')
output_high_contact = open('residue_lig_high_contact.txt', 'w')

#Compute % contact
per_contact = np.zeros(len(prot_res))
for i in range(len(prot_res)):
    dist_i = dist[i][:]
    len(dist_i)
    
    contact = 0 #conter for protein--ligand contact
    for j in dist_i:
        if j < 0.4:
            contact += 1
    per_contact[i] = 100 * (contact/len(dist_i))
    output.write(str(prot_res[i]) + ': ' + str(per_contact[i]))
    if per_contact[i] > 75:
        output_high_contact.write(prot_residue[i] + '\n')

