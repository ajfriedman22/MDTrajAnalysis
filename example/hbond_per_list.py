#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Percentage of time each individual h-bond is formed')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-p', required=True, help= 'Text file containing file path for h-bonds of interest (one file per line)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
File_path = args.p
miss_res = args.m

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 
import data_process 

sys.path.insert(1, prefix + '/protein_analysis/')
import hbond_analysis

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj)#Limit to uncorrelated frames
traj_prot = traj_uncorr.atom_slice(traj_uncorr.topology.select('protein')) #Select only atoms in the protein
top = traj_prot.topology
del traj; del traj_uncorr

#Set protein offset based on missing residues
offset = 1 + miss_res

#Make array for bond names from input file
name_bonds = []
input_paths = open(File_path, 'r').readlines()
for i in input_paths:
    file_path = i.strip()
    input_bonds = open(file_path, 'r').readlines()
    for n in input_bonds:
        n_clean = n.strip()
        if n_clean not in name_bonds:
            name_bonds.append(n_clean)

#Output file for bond 
output = open('Hbond_per_single.txt', 'w') #percentage of time each bond occurs

#Seperate residue names and atom names for each bond
res1, name1, res2, name2 = [],[],[],[]
for i in name_bonds:
    line = data_process.split(i.strip())
    #Determine indicise of dashed
    ind = []
    for j in range(len(line)):
        if line[j] == '-':
            ind.append(j)
    res1.append(str(int(data_process.sep_num(data_process.convert(line[0:ind[0]]))) - offset))
    name1.append(data_process.convert(line[ind[0]+1:ind[1]]).strip())
    res2.append(str(int(data_process.sep_num(data_process.convert(line[ind[2]+1:ind[3]]))) - offset))
    name2.append(data_process.convert(line[ind[3]+1:]).strip())

#Declare array for bond correlations
num_bonds = len(name_bonds)
bond_single_frac = np.zeros((num_bonds))

#track bond indicies
#Determine the percent of time each bond combination is present
for i in range(num_bonds):
    #If h-bond is between residues present in trajectory
    if int(res1[i]) >= 0 and int(res2[i]) >= 0 and int(res1[i]) < (traj_prot.n_residues-1) and int(res2[i]) < (traj_prot.n_residues-1):
        donor, acceptor, H = hbond_analysis.deter_bond(top, res1, res2, name1, name2, i)

        #Determine hydrogen with minimum distance
        H_min, dist = hbond_analysis.deter_H(acceptor, H, traj_prot)

        #Determine angle b/w donor and acceptor
        bond_a = np.array([[donor[0], H_min, acceptor[0]]])

        #Compute distances and angles over time for both bonds
        angle = md.compute_angles(traj_prot, bond_a , periodic = False)
        
        #Determine the percent of time both bonds are formed
        count_s=0 #single bonds
        for k in range(len(dist)):
            if dist[k] <= 0.25 and angle[k][0] >= 2.094:
                count_s += 1
        output.write(name_bonds[i] + ' ' + str(100*count_s/len(dist)) + '\n')
        bond_single_frac[i] = count_s/len(dist)
    else:
        output.write(name_bonds[i] + ' 0\n')

print('Hbond Percentages Calculated')
