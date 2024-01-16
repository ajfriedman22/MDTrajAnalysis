def com_rmsd(ref, traj, lig):
    import mdtraj as md
    import math
    import numpy as np

    #Align trajectory to reference
    traj_align = traj.superpose(ref)

    #seperate ligand carbon atoms
    lig_only_ref = ref.atom_slice(ref.topology.select('resname ' + str(lig))) #reference
    lig_only_traj = traj_align.atom_slice(traj_align.topology.select('resname ' + str(lig))) #trajectory

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

def deter_multimodal(dihedrals, name, i):
    import sys
    import os.path
    current_directory = os.path.dirname(os.path.realpath(__file__))
    prefix = current_directory.rsplit('/',1)[0]
    sys.path.insert(1, prefix + '/Traj_process/')
    import data_process

    #Seperate dihedral angles
    dihe_dist = dihedrals[:,i]
    
    #Determine maxima for probability distribution
    maxima = data_process.compute_max(dihe_dist)

    #Determine data not in the main peak
    main_peak, other_peak = [], []
    for i in dihe_dist:
        if abs(i - maxima) < 40 or abs(i + 360 - maxima) < 40 or abs(i - 360 - maxima) < 40:
            main_peak.append(i)
        else:
            other_peak.append(i)
    all_maxima = [data_process.compute_max(main_peak)]

    #If greater than 5% outliers count as seperate peak
    while len(other_peak)/len(dihe_dist) > 0.15:
        maxima = data_process.compute_max(other_peak)
        new_dist = other_peak
        main_peak, other_peak = [], []
        for i in new_dist:
            if abs(i - maxima) < 40 or abs(i + 360 - maxima) < 40 or abs(i - 360 - maxima) < 40:
                main_peak.append(i)
            else:
                other_peak.append(i)
        if len(main_peak) > (0.10*len(dihe_dist)):
            all_maxima.append(data_process.compute_max(main_peak))
        else:
            break
    return all_maxima, dihe_dist

