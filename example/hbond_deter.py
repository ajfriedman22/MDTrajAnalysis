#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination all h-bonds present more than set percent of the trajectory')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-f', required=False, type=int, default = 0.6, help= 'Minimum frequency h-bond appears(0 to 1)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
freq_set = args.f

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_analysis/')
import hbond_analysis

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj)#Limit to uncorrelated frames
traj_ns = traj_uncorr.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
del traj; del traj_uncorr #Save space

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine list of H-bonds present in the trajectory for over 60% of the frames
hbonds = md.baker_hubbard(traj_ns, freq=freq_set, exclude_water=True, periodic=False)
label = lambda hbond : '%s -- %s' % (traj_ns.topology.atom(hbond[0]), traj_ns.topology.atom(hbond[2])) #Extract labels for h-bonds
np.savetxt('Hbonds_atoms.txt', hbonds) #Save all atom indicies of h-bonds to file for further analysis

#Write all h-bonds present for >60% of trajectory to file
file_object = open('Hbonds_name.txt', 'w') 

for hbond in hbonds:
    file_object.write(label(hbond)) #Maintain same order as atom indicies
    file_object.write('\n')
file_object.close() #close file

#Determine the exact percentage of time that each h-bond present for >60% of the trajectory is formed
per = hbond_analysis.bond_per(traj_ns, hbonds)
np.savetxt('Hbonds_per.txt', per)
print('Hbond Percentages Calculated and Output')
