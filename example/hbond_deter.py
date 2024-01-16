#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
import pandas as pd

def split_res(res_atom, res_name_ar, res_num_ar, atom_ar):
    import re
    [res, atom] = res_atom.split('-')
    [res_name, res_num, empty] = re.split('(\d+)', res)
    res_name_ar.append(res_name)
    res_num_ar.append(res_num)
    atom_ar.append(atom)

#Silence MDTraj Warnings
warnings.filterwarnings("ignore")

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination all h-bonds present more than set percent of the trajectory')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-f', required=False, type=float, default = 0.6, help= 'Minimum frequency h-bond appears(0 to 1)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
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
traj = load_data.mdtraj_load(File_traj, File_gro, True, True)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine list of H-bonds present in the trajectory for over 60% of the frames
hbonds = md.baker_hubbard(traj, freq=freq_set, exclude_water=True, periodic=False)
label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])) #Extract labels for h-bonds

#Save atom indices to DataFrame
atom_df = pd.DataFrame({'Donor': hbonds[:,0], 'Hydrogen': hbonds[:,1], 'Acceptor': hbonds[:,2]})
atom_df.to_csv('Hbonds_atoms.csv')

#Write all h-bonds present for >60% of trajectory to file
donor_res_name, donor_res_num, donor_atom_name, acceptor_res_name, acceptor_res_num, acceptor_atom_name = [],[],[],[],[],[]
for hbond in hbonds:
    name = label(hbond) #Maintain same order as atom indicies
    [donor, acceptor] = name.split(' -- ')
    split_res(donor, donor_res_name, donor_res_num, donor_atom_name)
    split_res(acceptor, acceptor_res_name, acceptor_res_num, acceptor_atom_name)

#Determine the exact percentage of time that each h-bond present for >60% of the trajectory is formed
per = hbond_analysis.bond_per(traj, hbonds)

#Save names and percentages to file
df = pd.DataFrame({'Donor Residue Name': donor_res_name, 'Donor Residue ID': donor_res_num, 'Donor Atom Name': donor_atom_name, 'Acceptor Residue Name': acceptor_res_name, 'Acceptor Residue ID': acceptor_res_num, 'Acceptor Atom Name': acceptor_atom_name, 'Occupancy %': per})
df.to_csv('Hbond_per.csv')
print('Hbond Percentages Calculated and Output')
print('-------------------------------------------------------------')
