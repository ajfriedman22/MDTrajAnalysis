#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.constants import k, Avogadro
import random
import matplotlib.pyplot as plt
import seaborn as sns

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
    if np.abs(array[idx]-value) <= np.abs(array[idx2] - value + 360) and np.abs(array[idx] - value) <= np.abs(array[idx3] - value - 360):
        return idx+1
    elif np.abs(array[idx2] - value + 360) <= np.abs(array[idx] - value) and np.abs(array[idx2] - value + 360) <= np.abs(array[idx3] - value - 360):
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

def clust_conf(traj, per, file_name):
    #Compute pairwise distances
    distances = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        distances[i] = md.rmsd(traj, traj, i, atom_indices=traj.topology.select('element != H'))

    #Perform Clustering
    reduced_distances = squareform(distances, checks=False)
    link = linkage(reduced_distances, method='single') #The hierarchical clustering encoded as a matrix
    frame_list = dendrogram(link, no_labels=False, count_sort='descendent')['leaves']
    frame_cat = dendrogram(link, no_labels=False, count_sort='descendent')['color_list']

    #Keep only one file per cluster
    frames_sep = [] #List of frames that are unique and will be processed
    cat = frame_cat[0]
    frames_indv = [0]
    for frame in range(1, len(frame_list)-1):
        if frame_cat[frame] == cat:
            frames_indv.append(frame)
        else:
            frames_sep.append(frames_indv)
            cat = frame_cat[frame]
            frames_indv = [frame]
    frames_sep.append(frames_indv)

    #Save file names which have unique clusters
    frames_unique = []
    per_unique = np.zeros(len(frames_sep))
    for i in range(len(frames_sep)):
        if len(frames_sep[i]) > 0:
            for f in frames_sep[i]:
                per_unique[i] += per[f]
            if per_unique[i] > 0:
                frames_unique.append(int(random.sample(frames_sep[i], 1)[0]))
        else:
            raise Exception(f'Blank Cluster {frames_sep}')
    traj_clust = traj.slice(frames_unique)
    traj_clust.save_pdb(f'{file_name}.pdb')
    return traj_clust, per_unique, frames_sep

def process_confs(traj, per, file_name):
    #Save PDB
    traj.save_pdb(f'{file_name}.pdb')
    #Compute relative conformer energy
    rel_ener = get_rel_ener(per)

    #compute radius of gyration
    rg = md.compute_rg(traj)
    df_clust = pd.DataFrame({'Conformer ID': np.linspace(1, len(per), num=len(per), dtype=int), 'Occupancy': per, 'Relative FE': rel_ener})
    #print(df_clust)
    df_non_zero = df_clust[df_clust['Occupancy'] > 0]
    df_non_zero['Radius of Gyration'] = rg
    df_non_zero.to_csv(f'{file_name}.csv')

    labels = []
    for i, per in enumerate(df_non_zero['Occupancy']):
        if per > 1.5:
            labels.append(df_non_zero['Conformer ID'].values[i])
        else:
            labels.append('')
    plt.figure()
    plt.pie(df_non_zero['Occupancy'], labels=labels)
    plt.title(name)
    plt.savefig(f'{file_name}_pie.png')
    plt.close()

    plt.figure()
    sns.barplot(df_non_zero, x='Conformer ID', y='Relative FE')
    plt.ylabel('Relative Free Energy (kJ/mol)')
    plt.savefig(f'{file_name}_FE.png')
    plt.close()

    plt.figure()
    sns.barplot(df_non_zero, x='Conformer ID', y='Radius of Gyration')
    plt.ylabel(r'Radius of Gyration($\AA$)')
    plt.savefig(f'{file_name}_rg.png')
    plt.close()

    plt.figure()
    sns.histplot(df_non_zero, x='Radius of Gyration')
    plt.xlabel(r'Radius of Gyration($\AA$)')
    plt.savefig(f'{file_name}_rg_hist.png') 
    plt.close()

def get_rel_ener(per_all):
    per_non_zero = per_all[per_all != 0]
    ref_per = np.min(per_non_zero)
    rel_ener = np.zeros(len(per_all))
    for p, per in enumerate(per_all):
        rel_ener[p] = (k/1000)*300*np.log(per/ref_per)*Avogadro
    return rel_ener

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Ligand Conformers')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-s', required=True, type = str, help= 'name res# name_atom1 name_atom2 name_atom3 name_atom4')
parser.add_argument('-n', required=False, default = 'Lig', type = str, help= 'Ligand Name')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
file_input = args.s
name = args.n

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro, True, True)

#Cluster all conformers
traj_clust, per_clust, group = clust_conf(traj, np.ones(traj.n_frames), f'{name}_clust')
process_confs(traj_clust, per_clust, f'{name}_clust')

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
frame_select = []
for t in range(traj.n_frames):
    for i, conf in enumerate(conformer):
        if (conf == dihe_peak_sampled[t,:]).all():
            if count[i] == 0:
                frame_select.append(t)
            count[i] += 1
            break
per = 100*(count/traj.n_frames)

#Print conformer angle combinations, percent ligand is in conformation, and frame in which the ligand is in that conformation
df = pd.DataFrame({'Conformer ID': np.linspace(1, len(per), num=len(per), dtype=int)})
for i in range(num_dihe):
    df[f'Max for d{i+1}'] = conformer[:,i]
df.to_csv(f'{name}_dihe_def.csv')
traj_confs = traj.slice(frame_select)
process_confs(traj_confs, per, f'{name}_dihe')

#Cluster conformers
traj_dihe_clust, per_dihe_clust, group = clust_conf(traj_confs, per, f'{name}_dihe_clust')
df = pd.DataFrame({'Conformer ID': np.linspace(1, len(per_dihe_clust), num=len(per_dihe_clust), dtype=int), 'Grouped Confs': group})
df.to_csv(f'{name}_dihe_clust_def.csv')
process_confs(traj_dihe_clust, per_dihe_clust, f'{name}_dihe_clust')
