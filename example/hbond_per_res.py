#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
import pandas as pd

class ProgressBar(object):
    DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
    FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'

    def __init__(self, total, width=40, fmt=DEFAULT, symbol='=',
                 output=sys.stderr):
        assert len(symbol) == 1

        self.total = total
        self.width = width
        self.symbol = symbol
        self.output = output
        self.fmt = re.sub(r'(?P<name>%\(.+?\))d',
            r'\g<name>%dd' % len(str(total)), fmt)

        self.current = 0

    def __call__(self):
        percent = self.current / float(self.total)
        size = int(self.width * percent)
        remaining = self.total - self.current
        bar = '[' + self.symbol * size + ' ' * (self.width - size) + ']'

        args = {
            'total': self.total,
            'bar': bar,
            'current': self.current,
            'percent': percent * 100,
            'remaining': remaining
        }
        print('\r' + self.fmt % args, file=self.output, end='')

    def done(self):
        self.current = self.total
        self()
        print('', file=self.output)

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Percentage of time h-bond is formed between two given residues')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-p', required=True, help= 'Text file containing file path for h-bonds of interest (one file per line)')
parser.add_argument('-n', required=False, type=str, default='interest', help= 'Name for output file')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
File_path = args.p
miss_res = args.m
name = args.n

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 
import data_process 

sys.path.insert(1, prefix + '/protein_analysis/')
import hbond_analysis

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro, True, True)

#Set protein offset based on missing residues
offset = 1 + miss_res

#Make array for bond names from input file
name_bonds, options_bond = [],[]
input_paths = open(File_path, 'r').readlines()
for i in input_paths:
    input_df = pd.read_csv(i.strip())
    for index, row in input_df.iterrows():
        bond = str(row['Donor Residue ID']) + str(row['Donor Residue Name']) + ' - ' + str(row['Acceptor Residue ID']) + str(row['Acceptor Residue Name'])
        if bond not in name_bonds:
            name_bonds.append(bond)
            options_bond.append([row['Donor Residue Name'],  row['Donor Residue ID'], [row['Donor Atom Name']], row['Acceptor Residue Name'], row['Acceptor Residue ID'], [row['Acceptor Atom Name']]])
        else:
            if len(options_bond) == 0:
                options_bond.append([row['Donor Residue Name'],  row['Donor Residue ID'], [row['Donor Atom Name']], row['Acceptor Residue Name'], row['Acceptor Residue ID'], [row['Acceptor Atom Name']]])
            else:
                for n in range(len(options_bond)):
                    [d_res, d_id, d_atom, a_res, a_id, a_atom] = options_bond[n]
                    if row['Donor Residue ID'] == d_id and row['Acceptor Residue ID'] == a_id:
                        if row['Donor Atom Name'] not in [d_atom] or row['Acceptor Atom Name'] not in [a_atom]:
                            d_atom = [d_atom] + [row['Donor Atom Name']]
                            a_atom = [a_atom] + [row['Acceptor Atom Name']]
                            options_bond[n] = [d_res, d_id, d_atom, a_res, a_id, a_atom]

#Determine the percent of time each bond combination is present
per = np.zeros(len(name_bonds))
for i in range(len(name_bonds)):
    donor_res = int(options_bond[i][1]) - offset
    acceptor_res = int(options_bond[i][4]) - offset
    donor_options = options_bond[i][2]
    acceptor_options = options_bond[i][5]
    count_bond = np.zeros(traj.n_frames)

    #If h-bond is between residues present in trajectory
    if int(donor_res) >= 0 and int(acceptor_res) >= 0 and int(donor_res) < (traj.n_residues-1) and int(acceptor_res) < (traj.n_residues-1):
        for x in range(len(donor_options)):
            donor, acceptor, H = hbond_analysis.deter_bond(traj.topology, donor_res, acceptor_res, donor_options[x], acceptor_options[x])

            #Determine hydrogen with minimum distance
            H_min, dist = hbond_analysis.deter_H(acceptor, H, traj)

            #Determine angle b/w donor and acceptor
            bond_a = np.array([[donor[0], H_min, acceptor[0]]])

            #Compute distances and angles over time for both bonds
            angle = md.compute_angles(traj, bond_a , periodic = False)
        
            #Determine the percent of time both bonds are formed
            for k in range(len(dist)):
                if dist[k] <= 0.25 and angle[k][0] >= 2.094 and count_bond[k] == 0:
                    count_bond[k] = 1
        per[i] = 100*sum(count_bond)/len(dist)
df = pd.DataFrame({'Bond Name': name_bonds, 'Occupancy %': per})
df.to_csv('Hbonds_per_' + name + '.csv')

print('Hbond Percentages Calculated')
