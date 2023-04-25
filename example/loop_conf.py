import mdtraj as md
import numpy as np
import pandas as pd
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, protein and ligand RMSD, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-s', required=True, help= 'Input File(txt) Format: name res1 res2 dist_threshold')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
miss_res = args.m
input_file = args.s
if input_file.split('.')[-1] != 'txt': #Add default file extension if not in input
    input_file = input_file + '.txt'

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)

#Load and format input
input_data = open(input_file, 'r').readlines()
loop_name = []
loop_threshold = np.zeros(len(input_data))
loop_res = np.zeros((len(input_data), 2))
for i in range(len(input_data)):
    line = input_data[i].split()
    loop_name.append(line[0])
    loop_threshold[i] = float(line[3])
    loop_res[i,:] = [float(line[1]), float(line[2])]

#Calculate distances
dist = md.compute_contacts(traj, loop_res, scheme='closest-Heavy')
print(dist)

#Determine percent open
per = np.zeros(len(loop_name))
df = pd.DataFrame()
for n in range(len(loop_name)):
    count = 0
    loop_orient = []
    for t in range(traj.n_frames):
        if dist[t][n] < loop_threshold[n]:
            count += 1
            loop_orient.append(1)
            loop_orient.append(0)
    per[n] = 100*(count/traj.n_frames)
    df[str(loop_name[n])] = loop_orient
df.to_csv('loop_orientation_time.csv') 

#Save to file
df = pd.DataFrame({'Loop Name': loop_name, 'Percent Open': per})
df.to_csv('loop_orientation.csv')

print('Loop Oriantation Analysis Complete')
