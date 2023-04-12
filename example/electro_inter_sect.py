#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination all h-bonds present more than set percent of the trajectory')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-f', required=False, type=float, default = 0.6, help= 'Minimum frequency h-bond appears(0 to 1)')
parser.add_argument('-sect', required=False, type=str, default = 'none', help= 'File name for sections of interest')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
    File_traj = File_traj + '.xtc'
File_gro = args.g
if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
    File_gro = File_gro + '.gro'
miss_res = args.m
freq_set = args.f
input_file = args.sect

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_analysis/')
import hbond_analysis

#Load Trajectory
traj = load_data.mdtraj_load(File_traj, File_gro)
traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj)#Limit to uncorrelated frames
traj_ns = traj_uncorr.remove_solvent() #Remove solvent from the trajectory leaving only protein (and ligand if applicable)
del traj; del traj_uncorr #Save space

#Set protein offset based on missing residues
offset = 1 + miss_res

#Determine list of H-bonds present in the trajectory for over 60% of the frames
hbonds = md.baker_hubbard(traj_ns, freq=freq_set, exclude_water=True, periodic=False)
label = lambda hbond : '%s -- %s' % (traj_ns.topology.atom(hbond[0]), traj_ns.topology.atom(hbond[2])) #Extract labels for h-bonds

print('Hbonds Determined')

#Load protein section of interest
if input_file != 'none':
    sections = open(input_file, 'r').readlines()
    res_interest = []
    for i in range(len(sections)):
        name, sect = load_data.read_sections(sections, i, miss_res, 1, 1)
        for n in sect:
            res_interest.append(int(n))
else:
    res_interest = np.linspace(0, traj_ns.n_residues(), num = traj_ns.n_residues())

#Open output files
output_hbond = open('Hbonds_sect.txt', 'w')

#Determine % bond formed and sort salt bridges
da_distances = md.compute_distances(traj_ns, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
da_angles = md.compute_angles(traj_ns, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
[num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate
for j in range(len(hbonds)):
    bond = label(hbonds[j])
    res1, atom_name1, name1, res2, atom_name2, name2 = hbond_analysis.hbond_read(bond, offset)
    #Calculate the percent formed only if in desired range
    if int(res1) in res_interest or int(res2) in res_interest:
        count = 0 #Initialize count
        for i in range(num_t): #Loop through all frames
            if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If distance between donor and acceptor is less than 2.5A and the angle is greater than 120 degrees or ~ 2.094 radians
                count +=1
        per = 100*count/num_t #Percentage of time each h-bond is present in trajectory
        
        output_hbond.write(bond + ': ' + str(per) + '\n')

