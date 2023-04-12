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
parser.add_argument('-lref', required=False, type=str, default = 'none', help='Reference structure for Ligand COM RMSD')

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
lig_ref = args.lref

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_inter/')
import process_traj

sys.path.insert(1, prefix + '/ligand_analysis/')
import lig_motion

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_ns = traj.atom_slice(traj.topology.select('backbone or resname ' + lig_name))
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj_ns)#Limit to uncorrelated frames
del traj; del traj_ns

#Determine indices of protein and ligand residues in topology
lig_res = load_data.lig_check(lig, miss_res, traj_uncorr, lig_name)

#Load residues which define bound state
sections = open(res_bind_file, 'r').readlines()
res_bind = []
for i in range(len(sections)):
    name, sect = load_data.read_sections(sections, i, miss_res, 1, 1)
    for n in sect:
        res_bind.append(n)

#Determine distance between all protein residues and ligand
res_pairs = list(product([lig_res], res_bind))
[dist, pairs] = md.compute_contacts(traj_uncorr, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)

#Determine the % bound and the time ligand becomes unbound
frames, pairs = np.shape(dist)
lig_bind = np.zeros(frames) #Kepp track of if the ligand is bound or unbound in each frame
n_unbound = 0
frame_unbind = 'none'
for t in range(frames):
    n_inter = 0 #Count the number of interactions with residues defining binding
    for i in range(pairs):
        if dist[t][i] < 0.5:
            n_inter += 1
    if n_inter > 1:
        lig_bind[t] = 1
        n_unbound = 0
    else:
        lig_bind[t] = 0
        n_unbound += 1
    if n_unbound > 10*(frames/tot_time) and frame_unbind == 'none':
        frame_unbind = t - n_unbound
#Adjust for immediate unbinding
if frame_unbind != 'none' and int(frame_unbind) < 1:
    frame_unbind = (5/tot_time) * frames
per_bound = 100*(sum(lig_bind)/frames)
if frame_unbind == 'none':
    t_unbind = tot_time
    frame_unbind = frames
else:
    t_unbind = (frame_unbind/frames) * tot_time

output = open('Ligand_bind.txt', 'w')
output.write('Ligand bound for ' + str(per_bound) + '\n')
output.write('Unbinding time of ' + str(t_unbind) + '\n')
print('Ligand binding analysis completed') 

#Compute Ligand COM RMSD
if lig_ref != 'none':
    #Load reference
    ref_ns = load_data.load_ref(lig_ref, 'backbone or resname ' + str(lig_name))
    ref_top = ref_ns.topology

    #Limit trajectory to frames with bound ligand
    traj_bound = traj_uncorr[:frame_unbind]
    traj_ns = traj_bound.atom_slice(traj_uncorr.topology.select('backbone or resname ' + str(lig_name)))
    top_ns = traj_ns.topology

    traj_ns_align = traj_ns.superpose(ref_ns)
    
    #Compute COM displacment
    displacment, lig_rmsd = lig_motion.com_rmsd(ref_ns, traj_ns_align, lig_name)

    #Output COM RMSD
    output = open('lig_com_rmsd_' + lig_name + '.txt', 'w')
    output.write(str(lig_rmsd))

    #Output displacment for bootstrapping
    np.savetxt('lig_com_dis_' + lig_name + '.txt', displacment)

    print('Ligand COM RMSD Completed')
