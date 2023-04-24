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
    torsion_name = []
    for i in range(len(input_ind)):
        line = input_ind[i].split()
        torsion_name.append(line[0])
        for j in range(4):
            torsion_ind[i,j] = traj.topology.select('resid ' + str(int(line[1])-offset) + ' and name ' + str(line[j+2]))
    return torsion_name, torsion_ind

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Ligand Conformers')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-s', required=True, type = str, help= 'name res# name_atom1 name_atom2 name_atom3 name_atom4')
parser.add_argument('-n', required=True, type = str, help= 'Name for dihedral file')

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
traj = load_data.mdtraj_load(File_traj, File_gro)

#Set protein offset based on missing residues
offset = 1

#Load atom indices for torisonal angles from file
torsion_name, torsion_ind = input_torsion(file_input, traj)

#Compute dihedral angles for ligand
dihedral = md.compute_dihedrals(traj, indices=torsion_ind)

#Convert to degree
dihedral = dihedral*(180/np.pi)

#Plot and print angle distribution
dihe_max, dihe_ind = [],[]
for i in range(len(torsion_name)):
    maxima, dihe_dist = lig_motion.deter_multimodal(dihedral, torsion_name, i)
    plot.plot_torsion(dihe_dist, torsion_name, i, maxima)
    #If multiple peaks add to dihe_max array
    if len(maxima) >= 2:
        dihe_max.append(maxima)
        dihe_ind.append(i)

#Determine the number of sampled ligand conformers from the torsion angles
num_conf = 2**len(dihe_ind)
dihe_name = []
conf_per = np.zeros(num_conf+1) #Percent of trajectory in each coformer
conf_frame = np.zeros(num_conf+1) #First frame the conformer is observed
frames, num_di = np.shape(dihedral)
for t in range(0, frames): #Loop fromgh frames from trajectory
    c = 0 #Keep track of conformer number 
    for i in range(len(dihe_ind)):
        dihe_name.append('D' + str(c+1))
        [max_main, max_second] = dihe_max[i]
        n = dihe_ind[i]
        if abs(max_main - dihedral[t,n]) < 30 or abs(max_main - dihedral[t,n] + 360) < 30 or abs(max_main - dihedral[t,n] - 360) < 30:
            c = c #c is uncahnges
        elif abs(max_second - dihedral[t,n]) < 30 or abs(max_second - dihedral[t,n] + 360) < 30 or abs(max_second - dihedral[t,n] - 360) < 30:
            c = int(c + (num_conf/(2**(i+1))))
        else:
            #Does not fit in any conformers
            c = -1
            break
    conf_per[c] += 1
    if conf_frame[c] == 0:
        conf_frame[c] = t+1
conf_per = conf_per*(100/frames)
print(len(dihe_name))
print(len(conf_per))
#Check that sum of frames ligand is in each conformations total to the number of frames
if sum(conf_per) < 99.5:
    raise Warning('Error! Not all frames accounted for!')
elif sum(conf_per) > 100:
    raise Exception('Error! Total percent greater than 100%')

#Print conformer angle combinations, percent ligand is in conformation, and frame in which the ligand is in that conformation
df = pd.DataFrame({'Dihedral Name': dihe_name, 'Occupancy(%)': conf_per, 'Frame': conf_frame})
df2 = pd.DataFrame({'Dihedral': dihe_ind, 'Max 1': dihe_max[:][0], 'Max 2': dihe_max[:][1]})
df.to_csv('conf_' + name + '_per.csv')
df2.to_csv('conf_' + name + '_id.csv')

print('Dihedral Analysis Complete')
