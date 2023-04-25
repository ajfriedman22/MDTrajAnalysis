#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
import pandas as pd

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
top = traj.topology

#Set protein offset based on missing residues
offset = 1 + miss_res

#Seperate section for DSSP calculation
input_file = open(sect, 'r').readlines()

for i in range(len(input_file)):
    [name, first, last] = input_file[i].split()
    first, last = int(first), int(last)
    traj_sect = traj.atom_slice(top.select(str(int(first)-offset) + ' <= resid and resid <= ' + str(int(last)-offset)))

    #Compute Phi and Psi angles for all residues in the a7 helix in all frames
    phi_ind, phi_angle = md.compute_phi(traj_sect, periodic = True, opt = True)
    psi_ind, psi_angle = md.compute_psi(traj_sect, periodic = True, opt = True)
    time, angles = np.shape(phi_angle) #determine the number of frames and the number of angles computed
    
    #Compute Secondary Structure of all Residues in the a7 helix using MDTraj function
    dssp_raw = md.compute_dssp(traj_sect, simplified=False) #Compute DSSP for all residues in the a7 helix for all trajectory frames
    frames, residues = np.shape(dssp_raw)

    #Get vector of residues computed
    res_name = np.linspace(first, last, num=residues, dtype=int)

    #Output to file
    df = pd.DataFrame()
    for i in range(residues):
        res_i_dssp = []
        for j in range(frames):
            if dssp_raw[j][i] == ' ': #Replace spaces with l if input selected
                res_i_dssp.append('l')
            else:
                res_i_dssp.append(dssp_raw[j][i])
        df[str(res_name[i])] = res_i_dssp
    df.to_csv('DSSP/DSSP_' + name + '.csv') 

    #Save psi and phi angles to files and overwrite if file is present
    df = pd.DataFrame()
    for i in range(angles):
        df['Phi Angle ' + str(i+1)] = phi_angle[:,i]
        df['Psi Angle ' + str(i+1)] = psi_angle[:,i]
    df.to_csv('DSSP/phi_psi_' + name + '.csv')
print('DSSP Calculated and File Written')
