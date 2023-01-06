#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Pi stacking interactions with ligand')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-l', required=False, type=str, default = 0, help= 'Ligand residue name')
parser.add_argument('-ind', required=False, default='default', help= 'File containing atom indices for ligand torsions(txt)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Missing N-terminal protein residues')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
lig = args.l
file_ind = args.ind
if file_ind.split('.')[-1] != 'txt' and file_ind != 'default': #Add default file extension if not in input
    file_ind = file_ind + '.txt'
miss_res = args.m

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/ligand_analysis/')
import lig_motion

sys.path.insert(1, prefix + '/display_data/')
import plot

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_ns = traj.remove_solvent()
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj_ns)#Limit to uncorrelated frames
del traj; del traj_ns

#Set protein offset based on missing residues
offset = 1

#Load atom indices for ring atoms
if file_ind != 'default':
    input_ind = open(file_ind, 'r').readlines()
    lig_sp2_ind = np.zeros(len(input_ind))
    for i in range(len(input_ind)):
        line = input_ind[i].split()
        lig_sp2_ind[i] = int(line[1])-offset
elif lig == 'AD' or lig == 'AB':
    #Name of SP2 Carbons for common ligands
    if lig == 'AD':
        sp2_name = ['C10', 'C11', 'C14', 'C15']
    else:
        sp2_name = ['C2', 'C3', 'C11', 'C12']
    #Determin indices of sp2 carbons
    lig_sp2_ind = np.zeros(len(sp2_name))
    for i in range(len(sp2_name)):
        atom_id = traj_uncorr.topology.select('resname ' + lig + ' and name ' + sp2_name[i])
        lig_sp2_ind[i] = int(atom_id[0])
else:
    print('Atom Indices for ligand SP2 Carbons Need Supplied')
#Check that atom indices are for ligand
lig_atom = traj_uncorr.topology.select('resname ' + lig)
min_lig = min(lig_atom)
max_lig = max(lig_atom)
if ((min_lig <= lig_sp2_ind) & (max_lig >= lig_sp2_ind)).all():
    print('All atoms in ligand')
else:
    print('Atoms outside range for ligand! Exiting immediately!')
    exit()

#Determine all possible protein ring indices
traj_ring = traj_uncorr.atom_slice(traj_uncorr.topology.select('resname PHE or resname TYR'))#Seperate residues capable of forming pi-stacking interactions
ring_ind = traj_uncorr.topology.select('resname PHE or resname TYR')
num_ring_res = traj_ring.n_residues#determine number of residues capable of forming pi-stacking interactions
del traj_ring

res_name = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CG'] #Names for atoms in ring
prot_ring_ind = np.zeros((num_ring_res, len(res_name)))
top = traj_uncorr.topology
resid = 0
ring_resid = []
for i in range(num_ring_res):
    #determine residue id
    for r in range(resid+1, traj_uncorr.n_residues):
        if top.select('resid ' + str(r))[0] in ring_ind:
            resid = r
            ring_resid.append(r + 1 + miss_res)
            break
    for j in range(len(res_name)):
        atom = top.select('resid ' + str(resid) + ' and name ' + str(res_name[j]))
        prot_ring_ind[i,j] = atom[0]

#Determine if the ligand had pi-stacking interactions with any residues in the protein in each frame
time = traj_uncorr.n_frames
pi_tot = np.zeros(time) #Total number of pi-stacking interactions in each frame
pi_res_per = np.zeros(num_ring_res)
for r in range(num_ring_res): #For a pi-stacking interaction to be present three distances must be less than 4A
    res_ring_ind = prot_ring_ind[r,:]
    pair = list(product(res_ring_ind, lig_sp2_ind))
    num_pair = len(res_ring_ind)*len(lig_sp2_ind)
    dist = md.compute_distances(traj_uncorr, pair)
    for t in range(time):
        sp2_contact = 0
        for j in range(num_pair):
            if dist[t, j] < 0.5:
                sp2_contact += 1
        if sp2_contact > 2:
            pi_tot[t] += 1 #Add to count of total pi-stacking
            pi_res_per[r] += 1
pi_res_per = 100 * pi_res_per/time

#Save to file
output = open('pi_stack_lig.txt', 'w')
for i in range(len(ring_resid)):
    output.write(str(ring_resid[i]) + ' ' + str(pi_res_per[i]) + '\n')

np.savetxt('pi_stack_lig_tot.txt', pi_tot)
