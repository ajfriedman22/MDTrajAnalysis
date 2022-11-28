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
    sys.path.insert(1, '../Traj_process/')
    import process_data

    #Seperate dihedral angles
    dihe_name = name[i]
    dihe_dist = dihedrals[:,i]
    
    #Determine maxima for probability distribution
    maxima = process_data.compute_max(dihe_dist)

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
        maxima_main = process_data.compute_max(main_peak)
        maxima_second = process_data.compute_max(second_peak)

        return [maxima_main, maxima_second], dihe_dist
    else:
        return [maxima], dihe_dist

