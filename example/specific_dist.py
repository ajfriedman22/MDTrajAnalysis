#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of distance between residue pairs')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-s', required=False, default = 'none', help= 'File containing sections of interest(res1 atom1 res2 atom2)')

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
if sect.split('.')[-1] != 'txt' and sect != 'none': #Add default file extension if not in input
    sect = sect + '.txt'

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
top = traj.topology

#Get atom indices
file_sect = open(sect, 'r').readlines()
atom_ind = np.zeros((len(file_sect), 2))
n=0
for line in file_sect:
    #Split input file line
    res1, atom1, res2, atom2 = line.strip('\n').split(' ')
    
    #Determine atom indices for the two points
    atom1_ind = top.select(f"resid {int(res1)-1-miss_res} and name {atom1}")
    atom2_ind = top.select(f"resid {int(res2)-1-miss_res} and name {atom2}")
    atom_ind[n,:] = [atom1_ind[0], atom2_ind[0]]
    n+=1

#Compute distances
dist = md.compute_distances(traj, atom_pairs=atom_ind)

#Seperate and save files
n=0
for line in open(sect, 'r').readlines():
    #Split input file line
    res1, atom1, res2, atom2 = line.strip('\n').split(' ')
    
    #Save file
    np.savetxt(f"{res1}_{atom1}_{res2}_{atom2}_dist.txt", np.array(dist[:,n]))
    n += 1
