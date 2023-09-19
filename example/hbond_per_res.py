#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
import pandas as pd
from tqdm import tqdm
 
#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Percentage of time h-bond is formed between two given residues')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-p', required=True, help= 'Input file (csv or txt)')
parser.add_argument('-f', required=False, default='mult', help='Does input file contain file multiple file paths (mult) or is just one file being used (single)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
File_path = args.p
miss_res = args.m

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_analysis/')
import hbond_analysis

#Make array for bond names from input file
output_df = pd.DataFrame(columns=['Bond Name', 'Acceptor Atom', 'Donor Atom'])
if args.f == 'mult':
    input_paths = open(File_path, 'r').readlines()
else:
    input_paths = [File_path]

for i in input_paths:
    #Determine if input is pre-processed CSV or txt file listing bonds
    file_name = i.strip()
    if file_name.split('.')[-1] == 'csv':
        input_df = pd.read_csv(file_name)
    elif file_name.split('.')[-1] == 'txt':
        input_df = load_data.hbond_txt_to_df(file_name)

    #Seperate distinct h-bonds
    for index, row in input_df.iterrows():
        #Bond Name format
        bond = f"{row['Donor Residue Name']}-{row['Donor Residue ID']}-{row['Acceptor Residue Name']}-{row['Acceptor Residue ID']}"
        acceptor_atom = row['Acceptor Atom Name'].strip()
        donor_atom = row['Donor Atom Name'].strip()

        #Determine if bond is in dataframe already
        row_add = output_df.loc[output_df['Bond Name'] == bond]
        if len(row_add.index) == 0: #Bond not in list
            output_df.loc[len(output_df.index)] = {'Bond Name': bond, 'Donor Atom': np.array(donor_atom), 'Acceptor Atom': np.array(acceptor_atom)}
        elif len(row_add.index) == 1: #Bond in list
            #Determine if donor and accetor atom pair is in df already
            prev_accept = row_add.loc[row_add.index[0]].at["Acceptor Atom"] #get accpetor atom list
            prev_donor = row_add.loc[row_add.index[0]].at["Donor Atom"] #get donor atom list
            
            if acceptor_atom in prev_accept and donor_atom in prev_donor: #Skip if exact duplicate
                continue
            else:
                accept = np.append(prev_accept, acceptor_atom)
                donor = np.append(prev_donor, donor_atom)
                index = output_df.index[output_df['Bond Name'] == bond]
                output_df.loc[index[0]].at['Acceptor Atom'] = np.array(accept)
                output_df.loc[index[0]].at['Donor Atom'] =  np.array(donor)
        else:
            raise Exception('Each bond should only appear once in final DF')
    del input_df

#Declare empty arrays
num_bonds = len(output_df.index)
per = np.zeros(num_bonds)
print(f'{num_bonds} Unique Bonds Found')

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro, True, True)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine the percent of time each bond combination is present
for i in tqdm(range(num_bonds)):
    #Determine residue number and atom name for donor and acceptor
    bond = output_df['Bond Name'].iloc[i]
    res1 = int(bond.split("-")[1]) - offset
    res2 = int(bond.split("-")[3]) - offset
    name1_options = output_df['Donor Atom'].iloc[i]
    name2_options = output_df['Acceptor Atom'].iloc[i]
    
    #If h-bond is between residues present in trajectory
    frames = traj.n_frames
    if int(res1) >= 0 and int(res2) >= 0 and int(res1) < (traj.n_residues) and int(res2) < (traj.n_residues):
        if len(np.atleast_1d(name1_options)) > 1:
            num_opt = len(name1_options)
        else:
            num_opt = 1
        dist_all = np.zeros((num_opt, frames))
        angle_all = np.zeros((num_opt, frames))
        
        for x in range(num_opt):
            if num_opt > 1:
                donor, acceptor, H = hbond_analysis.deter_bond(traj.topology, res1, res2, name1_options[x], name2_options[x])
            else:
                donor, acceptor, H = hbond_analysis.deter_bond(traj.topology, res1, res2, name1_options, name2_options)

            #Determine hydrogen with minimum distance
            H_min, dist_all[x,:] = hbond_analysis.deter_H(acceptor, H, traj)
            
            #Determine angle b/w donor and acceptor
            bond_a = np.array([[donor[0], H_min, acceptor[0]]])

            #Compute distances and angles over time for both bonds
            angle = md.compute_angles(traj, bond_a , periodic = False)
            angle_all[x,:] = angle[:,0]

        #Determine the percent of time both bonds are formed
        count = 0
        for t in range(frames):
            index = np.where(dist_all[:,t] <= 0.3)
            if len(index[0]) > 0 and any(angle_all[x,t] >= 2 for x in index[0]):
                count += 1

        per[i] = 100*count/frames
    else:
        per[i] = 0
output_df['Occupancy(%)'] = per
output_df.to_csv('Hbond_res.csv')

print('Hbond Percentages Calculated')
