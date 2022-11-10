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
parser.add_argument('-ref2', required=False, type=str, default = 'none', help='Reference structure for RMSD if two proteins present')
parser.add_argument('-aa', required=False, default=0, type=int, help= 'Number of atoms in receptor')
parser.add_argument('-l', required=False, default = 'none', type = str, help= 'Ligand name if ligand analysis should be performed')
parser.add_argument('-lref', required=False, type=str, default = 'none', help='Reference structure for Ligand COM RMSD')
parser.add_argument('-s', required=False, default = 'none', type = str, help= 'File name for residue ranges for computed RMSD (sect_name ref_res_initial ref_res_final traj_res_initial traj_res_final)')
parser.add_argument('-rn', required=False, default = 'self', help='Reference name for RMSD')
parser.add_argument('-ln', required=False, default = 'self', help='Reference name for Ligand RMSD')

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
ref = args.ref
ref2 = args.ref2
miss_ref = args.mref
aa_atom = args.aa
lig = args.l
lig_ref = args.lref
lig_name = args.ln
rmsd_sect = args.s
ref_name = args.rn


#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_inter/')
import process_traj

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
top = traj.topology
    
#Set protein offset based on missing residues
offset = 1 + miss_res

#Add label to files if the trajectory is not equilibrated
if eq_time == 0:
    name_add = ''
else:
    name_add = 'unequil'

#Compute full BB RMSD
if ref != 'none':
    #Load reference
    if aa_atom != 0:
        ref_bb = load_data.load_ref(ref, 'backbone and index <= ' + str(aa_atom))
        traj_bb = traj.atom_slice(top.select('backbone and index <= ' + str(aa_atom))) #Backbond atoms of protein 1 only
    else:
        ref_bb = load_data.load_ref(ref, 'backbone')
        traj_bb = traj.atom_slice(top.select('backbone')) #Backbond atoms of protein only
    #Calculate RMSF from reference structure
    rmsf_data = md.rmsf(traj_bb, ref_bb, parallel=True, precentered=False)

    #Calculate RMSD for full protein relative to reference structure
    rmsd_full_uncorr, t_full = process_traj.compute_rmsd(traj_bb, ref_bb)
    
    if aa_atom  != 0 and ref2 != 'none':
        ref2_bb = load_data.load_ref(ref, 'backbone and index > ' + str(aa_atom))
        traj2_bb = traj.atom_slice(top.select('backbone and index > ' + str(aa_atom))) #Backbond atoms of protein 1 only
        #Calculate RMSF from reference structure
        rmsf_data2 = md.rmsf(traj2_bb, ref2_bb, parallel=True, precentered=False)

        #Calculate RMSD for full protein relative to reference structure
        rmsd_full_uncorr2, t_full = process_traj.compute_rmsd(traj2_bb, ref2_bb)

    #Only save RMSD and RMSF values to file for equilibrated input trajectories
    if eq_time != 0 and ref_name == 'self':
        np.savetxt('rmsd_full_ref_' + str(ref_name) + '.txt', rmsd_full_uncorr) #save to text file
        np.savetxt('rmsf_ref_' + str(ref_name) + '.txt', rmsf_data) #save to text file
        if aa_atom != 0 and ref2 != 'none':
            np.savetxt('rmsd_full_ref_' + str(ref_name) + '_prot2.txt', rmsd_full_uncorr2) #save to text file
            np.savetxt('rmsf_ref_' + str(ref_name) + '_prot2.txt', rmsf_data2) #save to text file

    np.savetxt('uncorrelated_frames' + name_add + '.txt', t_full) #Save indices for uncorrelated frames to file

    print('Full BB RMSD Computed')

    #Delete unneeded arrays to save memory
    del rmsf_data; del rmsd_full_uncorr; del ref_bb
else:
    print('Full BB RMSD Skipped')

#Compute BB RMSD by sections
if rmsd_sect != 'none':
    #Load reference
    traj_bb = traj.atom_slice(top.select('backbone')) #Backbond atoms of protein 1 only
    ref_bb = load_data.load_ref(ref, 'backbone')
    ref_top = ref_bb.topology
    top_bb = traj_bb.topology

    #Loop through sections from input file
    sections = open(rmsd_sect, 'r').readlines()
    for i in sections:
        #Split line by spaces
        line = i.split(' ')
        if len(line) == 5:
            [name, ref_res_initial, ref_res_final, traj_res_initial, traj_res_final] = line
        elif len(line) == 3:
            [name, traj_res_initial, traj_res_final] = line
            ref_res_initial = traj_res_initial
            ref_res_final = traj_res_final
        else:
            print('Error in input file')
            exit()

        #Reference offset
        offset_ref = miss_ref + 1

        #Seperate sections
        ref_res = [int(ref_res_initial) - offset_ref, int(ref_res_final) - offset_ref]
        if len(line) == 3:
            traj_res = ref_res
        else:
            traj_res = [int(traj_res_initial) - offset, int(traj_res_final) - offset]

        ref_sect = ref_bb.atom_slice(ref_top.select(str(ref_res[0]) + ' <= resid and resid <= ' + str(ref_res[1]))) #Limit trajectory to the section of choice
        traj_sect = traj_bb.atom_slice(top_bb.select(str(traj_res[0]) + ' <= resid and resid <= ' + str(traj_res[1]))) #Limit trajectory to the section of choice
        
        #Compute RMSD for section of interest
        rmsd_sect_uncorr, t_sect = process_traj.compute_rmsd(traj_sect, ref_sect)
        print(len(rmsd_sect_uncorr))

        #Save RMSD to file
        np.savetxt('rmsd/rmsd_' + name + '_ref_' + str(ref_name) + '.txt', rmsd_sect_uncorr)
    print('Section BB RMSD Completed')
else:
    print('Section BB RMSD Skipped')

#Compute Ligand COM RMSD
if lig != 'none':
    #Load reference
    aa_atom_ref = aa_atom - 1
    ref_ns = load_data.load_ref(lig_ref, '(backbone and index <= ' + str(aa_atom) + ') or resname ' + lig)
    traj_ns = traj.atom_slice(top.select('(backbone and index <= ' + str(aa_atom) + ') or resname ' + lig))
    ref_top = ref_ns.topology
    top_ns = traj_ns.topology

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
    lig_only_traj = traj_ns_align.atom_slice(top_ns.select('resname ' + str(lig))) #trajectory

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
    output = open('lig_com_rmsd_' + lig_name + '.txt', 'w')
    output.write(str(lig_rmsd))
    
    #Output displacment for bootstrapping
    np.savetxt('lig_com_dis_' + lig_name + '.txt', displacment)

    print('Ligand COM RMSD Completed')
    #Delete unneeded arays
    del displacment

else:
    print('Ligand COM RMSD Skipped')

#Compute Ligand Heavy atom RMSD
if lig != 'none':
    #Load reference
    ref_bb = load_data.load_ref(lig_ref, 'backbone or resname ' + lig)
    ref_top = ref_bb.topology
    
    #Seperate Ligand heavy atoms
    ref_sect = ref_bb.atom_slice(ref_top.select("resname " + lig)) #Limit trajectory to the section of choice
    traj_sect = traj.atom_slice(top.select('resname ' + lig)) #Limit trajectory to the section of choice

    #Compute RMSD for section of interest
    rmsd_sect_uncorr, t_sect = process_traj.compute_rmsd(traj_sect, ref_sect)
    
    #Save RMSD to file
    np.savetxt('rmsd_lig_heavy_ref_' + lig_name + '.txt', rmsd_sect_uncorr)

    print('Ligand Heavy Atom RMSD Completed')

else:
    print('Ligand Heavy Atom RMSD Skipped')

