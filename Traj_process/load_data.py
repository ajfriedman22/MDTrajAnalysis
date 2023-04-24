#!/ usr / bin / env python

#Load MDTraj trajectory and process to output one full trajectory, one with no solvent, one with all protein residues, and one with protein bb atoms only
#Input: File_traj = GROMACS trajectory in XTC format, File_gro = GROMACS GRO file, a7_res = residues which make up the a7 helix
#Output: 
#traj_bb = MDTraj trajectory with only protein bb atoms
#traj_prot = MDTraj trajectory with all protein atoms and residues
#traj_ns = MDTraj trajectory with solvent molecules removed
#traj_a7 = MDTraj trajectory with only the atoms in the a7 helix
#miss_first = Returns True if the first residue of PTP1B is missing and all indices need to be subtracted by 1
def mdtraj_load(File_traj, File_gro, rm_solvent=True):
    #Import required packages
    import mdtraj as md
    if File_traj.split('.')[-1] != 'xtc': #Add file extension if not in input
        File_traj = File_traj + '.xtc'
    if File_gro.split('.')[-1] != 'gro': #Add default file extension if not in input
        File_gro = File_gro + '.gro'

    #Load trajectories
    traj = md.load(File_traj, top=File_gro)
    
    print('Trajectory Loaded')
    if rm_solvent == True:
        traj_ns = traj.remove_solvent()
    
        #Remove uncorrelated Frames 
        traj_uncorr = remove_uncorr('uncorrelated_frames.txt', traj_ns)
    else:
        #Remove uncorrelated Frames 
        traj_uncorr = remove_uncorr('uncorrelated_frames.txt', traj)
    return traj_uncorr


#Load float data from a 2-column file with a prefix
#Input:
#file_dir = Directory input file is stored in
#file_name = Name of input file
#covert = Whether a conversion from nm to A is necessary
#Output: x, y = vectors for each of the two columns in the input file
def col2_float_data(file_dir, file_name, convert):
    x,y = [],[]
    #Process file_name
    if file_name.split('.')[-1] != 'xvg': #Add file extension if not in input
        file_name = file_name + '.xvg'

    #Load data
    with open(file_dir + '/' + file_name) as f:
        for _ in range(18):
            next(f)
        for line in f:
            cols = line.split()
            x.append(float(cols[0]))
            if convert == True:
                y.append(float(cols[1])*10)
            else:
                y.append(float(cols[1]))
    return x, y

def load_ref(ref, selection):
    import mdtraj as md

    #Load reference PDB
    ref_pdb = md.load_pdb(ref)
    top_ref = ref_pdb.topology
    
    #Limit reference PDB to Protein bb atoms
    ref_sect = ref_pdb.atom_slice(top_ref.select(selection))
    
    return ref_sect

def read_sections(sections, i, miss_res, num_prot_res, num_sect):
    import numpy as np

    line = sections[i].split()
    if len(line) == 6:
        [name1, sect1_start, sect1_end, name2, sect2_start, sect2_end] = line
    elif len(line) == 3:
        [name1, sect1_start, sect1_end] = line
        name2 = 'rest'
        sect2_start = 1 
        sect2_end = num_prot_res
    else:
        print('Error in sections input file!')
        exit()

    #Compute distance between all residues in sect1 and sect2
    sect1_start = int(sect1_start)-1-miss_res
    sect1_end = int(sect1_end)-1-miss_res
    sect2_start = int(sect2_start)-1-miss_res
    sect2_end = int(sect2_end)-1-miss_res

    sect1 = np.linspace(sect1_start, sect1_end, num=sect1_end-sect1_start+1)
    sect2 = np.linspace(sect2_start, sect2_end, num=sect2_end-sect2_start+1)
    
    if num_sect == 1:
        return name1, sect1
    else:
        return name1, name2, sect1, sect2

def lig_check(lig, miss_res, traj_ns, name):
    lig_res = lig - 1 - miss_res
    traj_lig = traj_ns.topology.select('resid ' + str(lig_res) + ' and resname ' + name)
    if len(traj_lig) == 0:
        print('Error in Ligand residue ID! Exiting Immediately!')
        exit()
    return lig_res

def remove_uncorr(file_name, traj):
    import numpy as np
    import mdtraj as md
    import os.path

    #Limit trajectory to uncorrelated frames
    if os.path.exists(file_name):
        uncorr_ind_string = open(file_name, 'r').readlines()
        uncorr_ind = np.zeros(len(uncorr_ind_string), dtype=int)
        for j in range(len(uncorr_ind_string)):
            uncorr_ind[j] = int(j)
        traj_uncorr = traj.slice(uncorr_ind)
    else:
        traj_uncorr = traj
    
    return traj_uncorr
