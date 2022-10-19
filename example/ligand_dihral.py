#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product
import matplotlib.pyplot as plt

def plot_torsion(dihedrals, name, i):
    #Seperate dihedral angles
    dihe_name = name[i]
    dihe_dist = dihedrals[:,i]
    
    #Convert to degree
    dihe_dist = dihe_dist*(180/np.pi)

    #Histogram of the data
    n, bins, patches = plt.hist(dihe_dist, 30, density=True, facecolor='g', alpha=0.75)

    plt.xlabel('Torsional Angle(rad)')
    plt.ylabel('Probability')
    plt.xlim(-180, 180)
    plt.title('Histogram of Torsion Angle ' + dihe_name)
    plt.grid(True)
    plt.savefig('dihedrals/dihe_angle_' + dihe_name + '.png')
    plt.close()

    np.savetxt('dihedrals/dihe_angle_' + dihe_name + '.txt', dihe_dist)

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, protein and ligand RMSD, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-l', required=False, type=str, default = 0, help= 'Ligand residue name')
parser.add_argument('-ind', required=True, help= 'File containing atom indices for ligand torsions(txt)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
lig = args.l
file_ind = args.ind
if file_ind.split('.')[-1] != 'txt': #Add default file extension if not in input
    file_ind = file_ind + '.txt'

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)

#Set protein offset based on missing residues
offset = 1

#Load atom indices for torisonal angles from file
input_ind = open(file_ind, 'r').readlines()
torsion_ind = np.zeros((len(input_ind), 4))
torsion_name = []
for i in range(len(input_ind)):
    line = input_ind[i].split()
    torsion_name.append(line[0])
    torsion_ind[i,:] = [int(line[1])-offset, int(line[2])-offset, int(line[3])-offset, int(line[4])-offset]

#Check that atom indices are for ligand
lig_atom = traj.topology.select('resname ' + lig)
min_lig = min(lig_atom)
max_lig = max(lig_atom)
if ((min_lig <= torsion_ind) & (max_lig >= torsion_ind)).all():
    print('All atoms in ligand')
else:
    print('Atoms outside range for ligand! Exiting immediately!')
    exit()

#Compute dihedral angles for ligand
dihedral = md.compute_dihedrals(traj, indices=torsion_ind)

#Plot and print angle distribution
for i in range(len(input_ind)):
    plot_torsion(dihedral, torsion_name, i)
