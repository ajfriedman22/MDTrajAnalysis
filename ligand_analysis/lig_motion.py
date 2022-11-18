def com_rmsd(ref, traj, lig):
    import mdtraj as md
    import math
    import numpy as np

    #seperate ligand carbon atoms
    lig_only_ref = ref.atom_slice(ref.topology.select('resname ' + str(lig))) #reference
    lig_only_traj = traj.atom_slice(traj.topology.select('resname ' + str(lig))) #trajectory

    lig_only_ref_top = lig_only_ref.topology
    lig_only_traj_top = lig_only_traj.topology
        
    #Compute COM of ligand
    com = md.compute_center_of_mass(lig_only_traj)
    com_ref = md.compute_center_of_mass(lig_only_ref)
        
    #Compute displacment
    time, dim = np.shape(com)
    displacment = np.zeros(time)
    for j in range(time):
        displacment[j] = (com[j][0] - com_ref[0][0])**2 + (com[j][1] - com_ref[0][1])**2 + (com[j][2] - com_ref[0][2])**2

    lig_rmsd = math.sqrt(np.mean(displacment))
    
    return displacment, lig_rmsd
