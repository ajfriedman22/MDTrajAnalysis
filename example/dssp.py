#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP for MD Trajectory')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-s', required=False, default = 'none', type = str, help= 'File name for residue ranges for computed DSSP (name initial final)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
sect = args.s


#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 
import process_traj

sys.path.insert(1, prefix + '/protein_analysis/')
import prot_struct

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj)#Limit to uncorrelated frames
top = traj_uncorr.topology
del traj
#Set protein offset based on missing residues
offset = 1 + miss_res

#Seperate section for DSSP calculation
input_file = open(sect, 'r').readlines()
[name, first, last] = input_file[0].split()
traj_sect = traj_uncorr.atom_slice(top.select(str(int(first)-offset) + ' <= resid and resid <= ' + str(int(last)-offset)))

#Compute Phi and Psi angles for all residues in the a7 helix in all frames
phi_ind, phi_angle = md.compute_phi(traj_sect, periodic = True, opt = True)
psi_ind, psi_angle = md.compute_psi(traj_sect, periodic = True, opt = True)
time, angles = np.shape(phi_angle) #determine the number of frames and the number of angles computed
    
#Compute Secondary Structure of all Residues in the a7 helix using MDTraj function
dssp_list = md.compute_dssp(traj_sect, simplified=False) #Compute DSSP for all residues in the a7 helix for all trajectory frames
file_dssp = open('DSSP_'+ name + '.txt','w') #Create output file for DSSP and write over if file is present

#Replace blank spaces with 'l' for loop to avoid output confusion
dssp_res_mod = prot_struct.dssp_remove_space(dssp_list, 'l')

#Output DSSP to file
frame_uncorr, residue = dssp_uncorr.shape
for i in range(frame_uncorr): #Each row is a single frame
    for j in range(residue):#Each column is a residue in the a7 helix
        file_dssp.write(dssp_res_mod[i,j] + ' ')
    file_dssp.write('\n') #New line between each time frame
file_dssp.close() #close file
    
#Save psi and phi angles to files and overwrite if file is present
np.savetxt('phi_' + name + '.txt', phi_angle)
np.savetxt('psi_' + name + '.txt', psi_angle)

