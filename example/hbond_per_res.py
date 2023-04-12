#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Percentage of time h-bond is formed between two given residues')
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
name_bonds, options_bond = [],[]
input_paths = open(File_path, 'r').readlines()
for i in input_paths:
    file_path = i.strip()
    input_bonds = open(file_path, 'r').readlines()
    for n in input_bonds:
        n_clean = n.strip()
        line = n_clean.split(' ')
        if line[0] not in name_bonds:
            name_bonds.append(line[0])
            options_bond.append(line[1:])

#Output file for bond 
output = open('Hbond_per_res.txt', 'w') #percentage of time each bond occurs

#Declare array for bond correlations
num_bonds = len(name_bonds)

#track bond indicies
#Determine the percent of time each bond combination is present
for i in range(num_bonds):
    res1 = int(options_bond[i][0]) - miss_res
    res2 = int(options_bond[i][1]) - miss_res
    name1_options = [options_bond[i][2], options_bond[i][4], options_bond[i][6], options_bond[i][8]]
    name2_options = [options_bond[i][3], options_bond[i][5], options_bond[i][7], options_bond[i][9]]
    count_bond = np.zeros(traj_prot.n_frames)

    #If h-bond is between residues present in trajectory
    if int(res1) >= 0 and int(res2) >= 0 and int(res1) < (traj_prot.n_residues-1) and int(res2) < (traj_prot.n_residues-1):
        for x in range(len(name1_options)):
            if name1_options[x] != 'nan':
                donor, acceptor, H = hbond_analysis.deter_bond(top, res1, res2, name1_options[x], name2_options[x])

                #Determine hydrogen with minimum distance
                H_min, dist = hbond_analysis.deter_H(acceptor, H, traj_prot)

                #Determine angle b/w donor and acceptor
                bond_a = np.array([[donor[0], H_min, acceptor[0]]])

                #Compute distances and angles over time for both bonds
                angle = md.compute_angles(traj_prot, bond_a , periodic = False)
        
                #Determine the percent of time both bonds are formed
                for k in range(len(dist)):
                    if dist[k] <= 0.25 and angle[k][0] >= 2.094 and count_bond[k] == 0:
                        count_bond[k] = 1

        output.write(name_bonds[i] + ' ' + str(100*sum(count_bond)/len(dist)) + '\n')
    else:
        output.write(name_bonds[i] + ' 0\n')

print('Hbond Percentages Calculated')
