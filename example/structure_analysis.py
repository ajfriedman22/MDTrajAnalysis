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
parser.add_argument('-e', required=False, default=0, type=int, help= 'If input trajectory is equilibrated input equilibration time otherwise input 0(ns)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-ref', required=False, type=str, default = 'none', help='Reference structure for RMSD')
parser.add_argument('-mref', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues for the reference structure(default 0)')
parser.add_argument('-df', required=True, help= 'Directory github repo is in')
parser.add_argument('-l', required=False, default = 'none', type = str, help= 'Ligand name if ligand analysis should be performed')
parser.add_argument('-lref', required=False, type=str, default = 'none', help='Reference structure for Ligand COM RMSD')
parser.add_argument('-s', required=False, default = 'none', type = str, help= 'File name for residue ranges for computed RMSD (sect_name ref_res_initial ref_res_final traj_res_initial traj_res_final)')
parser.add_argument('-rn', required=False, default = 'self', help='Reference name for RMSD')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
eq_time = args.e
miss_res = args.m
ref_name = args.ref
miss_ref = args.mref
directory = args.df
lig = args.l
lig = args.lref
rmsd_sect = args.s
ref_name = args.rn

#Source custom functions
prefix = directory + '/MDTrajAnalysis/'
sys.path.insert(1, prefix + 'Traj_process/')
import load_data 

sys.path.insert(1, prefix + 'protein_inter/')
import process_traj

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Add label to files if the trajectory is not equilibrated
if eq_time == 0:
    name_add = ''
else:
    name_add = 'unequil'

#Compute full BB RMSD
if ref_name != 'none' and traj_bb.n_residues == ref_bb.n_residues:
    #Load reference
    ref_bb = load_data.load_ref(ref, 'backbone')

    #Calculate RMSF from reference structure
    rmsf_data = md.rmsf(traj_bb, ref_bb, parallel=True, precentered=False)

    #Calculate RMSD for full protein relative to reference structure
    rmsd_full_uncorr, t_full = process_traj.compute_rmsd(traj_bb, ref_bb, False)
        
    #Determine time in ns for each uncorrelated frame
    frame_per_ns = traj_bb.n_frames/traj_time
    uncorr_time = np.zeros(len(t_full))
    for i in range(len(t_full)):
        uncorr_time[i] = t_full[i]/frame_per_ns + eq_time

    #Only save RMSD and RMSF values to file for equilibrated input trajectories
    if eq_time != 0 and ref_name == 'self':
        np.savetxt('rmsd_full_ref_' + str(ref_name) + '.txt', rmsd_full_uncorr) #save to text file
        np.savetxt('rmsf_ref_' + str(ref_name) + '.txt', rmsf_data) #save to text file
    np.savetxt('uncorrelated_frames' + name_add + '.txt', t_full) #Save indices for uncorrelated frames to file
    np.savetxt('uncorrelated_time' + name_add + '.txt', uncorr_time) #Save time for uncorrelated frames to file

    print('Full BB RMSD Computed')

    #Delete unneeded arrays to save memory
    del rmsf_data; del rmsd_full_uncorr
else:
    print('Full BB RMSD Skipped')

#Compute BB RMSD by sections
if rmsd_sect != 'none':
    #Load reference
    ref_bb = load_data.load_ref(ref, 'backbone')
    ref_top = ref_bb.topology

    #Loop through sections from input file
    sections = open(rmsd_sect, 'r').readlines()
    for i in sections:
        #Split line by spaces
        [name, ref_res_initial, ref_res_final, traj_res_initial, traj_res_final] = i.split(' ')
        
        #Reference offset
        offset_ref = miss_ref + 1

        #Seperate sections
        ref_sect = ref_bb.atom_slice(ref_top.select(str(ref_res_initial - offset_ref) + ' >= resid and resid <= ' + str(ref_res_final - offset_ref))) #Limit trajectory to the section of choice
        traj_sect = traj_bb.atom_slice(traj_top.select(str(traj_res_initial - offset) + ' >= resid and resid <= ' + str(traj_res_final - offset))) #Limit trajectory to the section of choice

        #Compute RMSD for section of interest
        rmsd_sect_uncorr, t_sect = process_traj.compute_rmsd(traj_sect, ref_sect)
    
        #Save RMSD to file
        np.savetxt('rmsd_' + name + '_ref_' + str(ref_name) + '.txt', rmsd_sect_uncorr)
    print('Section BB RMSD Completed')
else:
    print('Section BB RMSD Skipped')

#Compute Ligand COM RMSD
if lig != 'none':
    #Load reference
    ref_ns = load_data.load_ref(ref, 'backbone or resname ' + lig)
    traj_ns = traj.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)

    #Limit trajectory to uncorrelated frames
    if os.path.exists('uncorrelated_frames.txt'):
        uncorr_ind_string = open('uncorrelated_frames.txt', 'r').readlines()
        uncorr_ind = np.zeros(len(uncorr_ind_string), dtype=int)
        for j in range(len(uncorr_ind_string)):
            uncorr_ind[j] = int(j)
        traj_uncorr = traj_ns.slice(uncorr_ind)
    else:
        traj_uncorr = traj_ns

    traj_ns_align = traj_uncorr.superpose(ref_ns)

    #seperate ligand carbon atoms
    lig_only_ref = ref_ns.atom_slice(ref_top.select('resname ' + str(lig))) #reference
    lig_only_traj = traj_ns_align.atom_slice(top.select('resname ' + str(lig))) #trajectory

    lig_only_ref_top = lig_only_ref.topology
    lig_only_traj_top = lig_only_traj.topology
        
    #Compute COM of ligand
    com = md.compute_center_of_mass(lig_only_traj)
    com_ref = md.compute_center_of_mass(lig_only_ref)
        
    #Compute displacment
    time, dim = np.shape(com)
    displacment = np.zeros(time)
    for j in range(time):
        displacment[j] = (com[j][0] - com_ref[0][0])**2 + (com[j][1] - com_ref[0][1])**2 + (com[j][2] - com_ref[0][2])**2

    lig_rmsd = math.sqrt(np.mean(displacment))
    
    #Output COM RMSD
    output = open('lig_com_rmsd.txt', 'w')
    output.write(lig_rmsd)
    
    #Output displacment for bootstrapping
    np.savetxt('lig_com_dis.txt', displacment)

    print('Ligand COM RMSD Completed')
    #Delete unneeded arays
    del displacment

else:
    print('Ligand COM RMSD Skipped')
