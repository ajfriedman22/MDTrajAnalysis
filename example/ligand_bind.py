#!/ usr / bin / env python
import math
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
parser.add_argument('-l', required=True, type = int, help= 'Ligand residue number')
parser.add_argument('-ln', required=False, default = 'LIG', type = str, help= 'Ligand name in GRO file')
parser.add_argument('-bind', required=True, help= 'File containing residues refingin bound state')
parser.add_argument('-n', required=True, type = int, help= 'Total time of trajectory(ns)')

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
lig_name = args.ln
res_bind_file = args.bind
if res_bind_file.split('.')[-1] != 'txt': #Add default file extension if not in input
    res_bind_file = res_bind_file + '.txt'
tot_time = args.n

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj)#Limit to uncorrelated frames
traj_ns = traj_uncorr.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
top = traj_ns.topology
del traj; del traj_uncorr
print('Trajectory Loaded')

#Determine indices of protein and ligand residues in topology
lig_res = load_data.lig_check(lig, miss_res, traj_ns, lig_name)

#Load residues which define bound state
sections = open(res_bind_file, 'r').readlines()
res_bind = []
for i in range(len(sections)):
    name, sect = load_data.read_sections(sections, i, miss_res, 1, 1)
    for n in sect:
        res_bind.append(n)

#Determine distance between all protein residues and ligand
res_pairs = list(product([lig_res], res_bind))
[dist, pairs] = md.compute_contacts(traj_ns, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)

#Determine the % bound and the time ligand becomes unbound
frames, pairs = np.shape(dist)
lig_bind = np.zeros(frames) #Kepp track of if the ligand is bound or unbound in each frame
n_unbound = 0
frame_unbind = 'none'
for t in range(frames):
    n_inter = 0 #Count the number of interactions with residues defining binding
    for i in range(pairs):
        if dist[t][i] < 0.4:
            n_inter += 1
    if n_inter > 2:
        lig_bind[t] = 1
        n_unbound = 0
    else:
        lig_bind[t] = 0
        n_unbound += 1
    if n_unbound > 10*(frames/tot_time) and frame_unbind == 'none':
        frame_unbind = t - n_unbound
per_bound = 100*(sum(lig_bind)/frames)
if frame_unbind == 'none':
    t_unbind = tot_time
else:
    t_unbind = (frame_unbind/frames) * tot_time

output = open('Ligand_bind.txt', 'w')
output.write('Ligand bound for ' + str(per_bound) + '\n')
output.write('Unbinding time of ' + str(t_unbind) + '\n')
print('Ligand binding analysis completed') 
