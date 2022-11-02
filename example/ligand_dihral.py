#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy import optimize
from scipy.signal import argrelextrema

def compute_max(data):
    kde = gaussian_kde(data)
    opt = optimize.minimize_scalar(lambda x: -kde(x))
    samples = np.linspace(min(data), max(data), 50)
    probs = kde.evaluate(samples)
    maxima_index = probs.argmax()
    maxima = samples[maxima_index]
    
    return maxima

def plot_torsion(dihe_dist, name, i, maxima):
    #Seperate dihedral angles
    dihe_name = name[i]
    
    #Histogram of the data
    n, bins, patches = plt.hist(dihe_dist, 30, density=True, facecolor='g', alpha=0.75)
    #Inidcate Maxima
    plt.axvline(x = maxima[0], color = 'k')
    if len(maxima) == 2:
        plt.axvline(x = maxima[1], color = 'b')

    plt.xlabel('Torsional Angle(rad)')
    plt.ylabel('Probability')
    plt.xlim(-180, 180)
    plt.title('Histogram of Torsion Angle ' + dihe_name)
    plt.grid(True)
    plt.savefig('dihedrals/dihe_angle_' + dihe_name + '.png')
    plt.close()
    
def deter_multimodal(dihedrals, name, i):
    #Seperate dihedral angles
    dihe_name = name[i]
    dihe_dist = dihedrals[:,i]
    
    #Determine maxima for probability distribution
    maxima = compute_max(dihe_dist)

    #Determine data not in the main peak
    main_peak = []
    second_peak = []
    for i in dihe_dist:
        diff = i - maxima
        if abs(i - maxima) < 30 or abs(i + 360 - maxima) < 30 or abs(i - 360 - maxima) < 30:
            main_peak.append(i)
        else:
            second_peak.append(i)
        
    #If greater than 3 outliers count as seperate peak
    if len(second_peak) > (0.05*len(dihe_dist)):
        #Determine the maxima of the two peaks individually 
        maxima_main = compute_max(main_peak)
        maxima_second = compute_max(second_peak)

        return [maxima_main, maxima_second], dihe_dist
    else:
        return [maxima], dihe_dist

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

#Convert to degree
dihedral = dihedral*(180/np.pi)

#Plot and print angle distribution
dihe_max, dihe_ind = [],[]
for i in range(len(input_ind)):
    maxima, dihe_dist = deter_multimodal(dihedral, torsion_name, i)
    plot_torsion(dihe_dist, torsion_name, i, maxima)
    #If multiple peaks add to dihe_max array
    if len(maxima) >= 2:
        dihe_max.append(maxima)
        dihe_ind.append(i)

#Determine the number of sampled ligand conformers from the torsion angles
num_conf = 2**len(dihe_ind)
conf_per = np.zeros(num_conf+1) #Percent of trajectory in each coformer
conf_frame = np.zeros(num_conf+1) #First frame the conformer is observed
frames, num_di = np.shape(dihedral)
for t in range(0, frames): #Loop fromgh frames from trajectory
    c = 0 #Keep track of conformer number 
    for i in range(len(dihe_ind)):
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

#Check that sum of frames ligand is in each conformations total to the number of frames
if sum(conf_per) < 99.5:
    print('Error! Not all frames accounted for!')
elif sum(conf_per) > 100:
    print('Error! Total percent greater than 100%')

#Print conformer angle combinations, percent ligand is in conformation, and frame in which the ligand is in that conformation
output_per = open('ligand_conf_per.txt', 'w')
for i in range(len(conf_per)):
    if conf_per[i] > 0:
        output_per.write('Conformation ' + str(i+1) + ': ' + str(conf_per[i]) + '%\n')
        output_per.write('frame ' + str(conf_frame[i]) + '\n')
    
output_conf = open('ligand_conf_id.txt', 'w')
for i in range(len(dihe_ind)):
    output_conf.write('Diherdal ' + str(dihe_ind[i]) + ': ' + str(dihe_max[i][0]) + ' or ' + str(dihe_max[i][1]) + '\n')
