#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

def hbond_deter(traj, file_name, freq):
    hbonds = md.baker_hubbard(traj, freq=freq, exclude_water=True, periodic=False)
    
    #Determine % hbond formed
    per = [] #Declare empty array for percentage of time h-bond is formed
    da_distances = md.compute_distances(traj, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
    da_angles = md.compute_angles(traj, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
    [num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate
    for j in range(num_h): #Loop through all h-bonds
        count = 0 #Initialize count
        for i in range(num_t): #Loop through all frames
            if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If distance between donor and acceptor is less than 2.5A and the angle is greater than 120 degrees or ~ 2.094 radians
                count +=1
        per.append(100*count/num_t) #Percentage of time each h-bond is present in trajectory

    label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2])) #Extract labels for h-bonds

    #Write all h-bonds present for >60% of trajectory to file
    file_object = open(file_name, 'w') 
    
    i=0
    for hbond in hbonds:
        file_object.write(label(hbond) + ': ' + str(per[i]) + '\n') #Maintain same order as atom indicies
        i+=1
    file_object.close() #close file

    return hbonds

def hbond_per(traj, file_name, hbonds):
    per = [] #Declare empty array for percentage of time h-bond is formed
    da_distances = md.compute_distances(traj, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
    da_angles = md.compute_angles(traj, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
    [num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate
    for j in range(num_h): #Loop through all h-bonds
        count = 0 #Initialize count
        for i in range(num_t): #Loop through all frames
            if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If distance between donor and acceptor is less than 2.5A and the angle is greater than 120 degrees or ~ 2.094 radians
                count +=1
        per.append(100*count/num_t) #Percentage of time each h-bond is present in trajectory
    np.savetxt(file_name, per)


#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of DSSP, H-bonds, Ligand Contacts, protein and ligand RMSD, Helical interactions and PCA for GROMACS Trajectory of PTP1B')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-l', required=False, type=int, default = 0, help= 'Ligand residue ID')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
lig = args.l

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
top = traj.topology

#Set protein offset based on missing residues
offset = 1 + miss_res

#Limit trajectory to ligand and protein residues of interest
#traj_lig = traj.atom_slice(top.select('resname LIG or water')) 
traj_sect = traj.atom_slice(top.select('resname LIG or (80 <= resid and resid <= 98) or (55 <= resid and resid <= 70) or (110 <= resid and resid >= 125) or (145 <= resid and resid <= 170)'))

#Determine all h-bonds present in trajectory b/w ligand and protein
#hbond_deter(traj_lig, 'Hbonds_lig_water.txt', 0.7)

#Determin all hbonds present b/w ligand and water
hbond_deter(traj_sect, 'Hbonds_lig_prot.txt', 0.1)



