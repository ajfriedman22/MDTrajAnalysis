#!/ usr / bin / env python

#Load MDTraj trajectory and process to output one full trajectory, one with no solvent, one with all protein residues, and one with protein bb atoms only
#Input: File_traj = GROMACS trajectory in XTC format, File_gro = GROMACS GRO file, a7_res = residues which make up the a7 helix
#Output: 
#traj_bb = MDTraj trajectory with only protein bb atoms
#traj_prot = MDTraj trajectory with all protein atoms and residues
#traj_ns = MDTraj trajectory with solvent molecules removed
#traj_a7 = MDTraj trajectory with only the atoms in the a7 helix
#miss_first = Returns True if the first residue of PTP1B is missing and all indices need to be subtracted by 1
def mdtraj_load(File_traj, File_gro):
    #Import required packages
    import mdtraj as md

    #Load trajectories
    traj = md.load(File_traj, top=File_gro)
    
    print('Trajectory Loaded')
    return traj

#Load float data from a 2-column file with a prefix
#Input:
#file_dir = Directory input file is stored in
#file_name = Name of input file
#covert = Whether a conversion from nm to A is necessary
#Output: x, y = vectors for each of the two columns in the input file
def col2_float_data(file_dir, file_name, convert):
    x,y = [],[]
    #Load data
    with open(file_dir + '/' + file_name + '.xvg') as f:
        for _ in range(18):
            next(f)
        for line in f:
            cols = line.split()
            x.append(float(cols[0]))
            if convert == True:
                y.append(float(cols[1])*10)
            else:
                y.append(float(cols[1]))
    return x,y

