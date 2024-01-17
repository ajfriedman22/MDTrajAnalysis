#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product
import pandas as pd

def input_torsion(file_input, traj):
    input_ind = open(file_input, 'r').readlines()
    torsion_ind = np.zeros((len(input_ind), 4))
    torsion_name, max_values, peak_options = [], [], []
    for i in range(len(input_ind)):
        line = input_ind[i].split()
        torsion_name.append(line[0])
        for j in range(4):
            torsion_ind[i,j] = traj.topology.select('resid ' + str(int(line[1])-offset) + ' and name ' + str(line[j+2]))
        max_values.append(line[6:])
        opt = np.linspace(1, len(line[6:]), num=len(line[6:]), dtype=int)
        peak_options.append(list(opt))
    return torsion_name, torsion_ind, max_values, peak_options

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    idx2 = (np.abs(array - value + 360)).argmin()
    idx3 = (np.abs(array - value - 360)).argmin()
    if array[idx] <= array[idx2] and array[idx] <= array[idx3]:
        return idx+1
    elif array[idx2] <= array[idx] and array[idx2] <= array[idx3]:
        return idx2+1
    else:
        return idx3+1
def get_conformations(peak_options):
    num_combos = len(peak_options[0])
    for i in range(1, len(peak_options)):
        num_combos = num_combos*len(peak_options[i])
    combos = np.zeros((num_combos, len(peak_options)))
    c = 0
    if len(peak_options) == 6:
        for i in peak_options[0]:
            for j in peak_options[1]:
                for k in peak_options[2]:
                    for l in peak_options[3]:
                        for m in peak_options[4]:
                            for n in peak_options[5]:
                                combos[c,:] = [i, j, k, l, m, n]
                                c+=1
    else:
        for i in peak_options[0]:
            for j in peak_options[1]:
                for k in peak_options[2]:
                    for l in peak_options[3]:
                        combos[c,:] = [i, j, k, l]
                        c +=1
    return combos
#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Ligand Conformers')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-s', required=True, type = str, help= 'name res# name_atom1 name_atom2 name_atom3 name_atom4')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
file_input = args.s

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
traj = load_data.mdtraj_load(File_traj, File_gro, True, True)

#Set protein offset based on missing residues
offset = 1

#Load atom indices for torisonal angles from file
torsion_name, torsion_ind, max_values, peak_options = input_torsion(file_input, traj)

#Compute dihedral angles for ligand
dihedral = md.compute_dihedrals(traj, indices=torsion_ind)

#Convert to degree
dihedral = dihedral*(180/np.pi)

#Determine which dihderal peak is being sampled per frame
num_dihe = len(torsion_name)
dihe_peak_sampled = np.zeros((traj.n_frames, num_dihe))
for t in range(traj.n_frames):
    for i in range(num_dihe):
        max_value_i = np.array(max_values[i], dtype=float)
        value = dihedral[t,i]
        dihe_peak_sampled[t,i] = find_nearest(max_value_i, value)
print('Peaks Found')

#Name conformations
conformer = get_conformations(peak_options)

#Classify dihedrals into conformations
count = np.zeros(len(conformer))
for t in range(traj.n_frames):
    for i, conf in enumerate(conformer):
        if (conf == dihe_peak_sampled[t,:]).all():
            count[i] += 1
            break
per = 100*(count/traj.n_frames)

#Print conformer angle combinations, percent ligand is in conformation, and frame in which the ligand is in that conformation
df = pd.DataFrame({'Conformer ID': np.linspace(1, len(per), num=len(per), dtype=int), 'Occupancy': per})
for i in range(1, num_dihe):
    df[f'Max for d{i}'] = conformer[:,int(i-1)]
df_non_zero = df[df['Occupancy'] > 0]
df_non_zero.to_csv('conf_per.csv')

print('Dihedral Analysis Complete')
print('---------------------------------------------------------------')
