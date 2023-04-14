#!/ usr / bin / env python
import math
import seaborn as sns
import pandas as pd
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
parser.add_argument('-ref', required=False, type=str, default = 'none', help='Reference structure for RMSD')
parser.add_argument('-mref', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues for the reference structure(default 0)')
parser.add_argument('-input', required=True, type=str, help= 'Text file containing residues and atom names for measuring helix twisting')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
ref = args.ref
miss_ref = args.mref
file_input = args.input

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_inter/')
import process_traj

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj)
traj_prot = traj_uncorr.atom_slice(traj_uncorr.topology.select('protein')) #Backbond atoms of protein 1 only
del traj; del traj_uncorr

#Set protein offset based on missing residues
offset = 1 + miss_res
ref_offset = 1 + miss_ref

#Load Reference
ref_prot = load_data.load_ref(ref, 'protein')

#Load input file and determine atom numbers
atom_txt = open(file_input, 'r').readlines()
atom_ref, atom_traj, atom_label, CA_ref = [],[],[],[]
for line in atom_txt:
    atom_num, atom_name = line.split(' ')
    atom_ref.append(ref_prot.topology.select('resi ' + str(int(atom_num) - ref_offset) + ' and name ' + str(atom_name)))
    CA_ref.append(ref_prot.topology.select('resi ' + str(int(atom_num) - ref_offset) + ' and name CA'))
    atom_traj.append(traj_prot.topology.select('resi ' + str(int(atom_num) - offset) + ' and name ' + str(atom_name)))

    atom_label.append(atom_num)

#Align reference to trajectory
traj_align = traj_prot.superpose(ref_prot, atom_indices=traj_prot.topology.select('resi > ' + str(186-offset) + ' and resi < ' + str(200-offset) + ' and backbone'), ref_atom_indices=ref_prot.topology.select('resi > ' + str(186-ref_offset) + ' and resi < ' + str(200-ref_offset) + ' and backbone'))

#Determine distance from reference for each atom pairs
res_dist = np.zeros((len(atom_label), traj_prot.n_frames))
iso_length = np.zeros(len(atom_label))
for n in range(len(atom_label)):
    #Determine corrdinates for reference points
    p1 = ref_prot.xyz[0, atom_ref[n], :]
    p2 = traj_align.xyz[:, atom_traj[n], :]
    p2 = np.squeeze(p2)

    #Determine distances b/w coordinates
    res_dist[n,:] = np.sqrt(np.sum((p1-p2)**2, axis=1))

    #Determine Length of residue to helix center for each residue of interest
    p1_CA = ref_prot.xyz[0, CA_ref[n],:]
    iso_length[n] = np.sqrt(np.sum((p1-p1_CA)**2, axis=1)) + (.54/2) #Triangle vertex in center of helix

print('Residue Distances Computed')

#Determine angle from residue distances
res_angle = np.zeros((len(atom_label), traj_prot.n_frames))
for n in range(len(atom_label)):
    res_angle[n,:] = np.arccos((2*iso_length[n]**2 - res_dist[n,:]**2)/(2*iso_length[n]**2))

    np.savetxt(atom_label[n] + '_angle.txt', res_angle)
