def traj_sect(traj, prot_res, offset):
    import mdtraj as md
    res = prot_res - offset
    top = traj.topology
    traj_sect = top.select('resid ' + str(res)) #Select only atoms in the given section
    return traj_sect

def water_neighbor(traj, prot_res, offset):
    import mdtraj as md
    import numpy as np

    #Seperate protein residues of interest
    prot_ind = traj_sect(traj, prot_res, offset)

    #Compute neighboring atoms for all residues of interest
    neighbors = md.compute_neighbors(traj, 0.5, prot_ind, haystack_indices=None, periodic=True)
    
    #Get list of all unique neighboring water atoms
    top = traj.topology
    traj_solv = top.select('water')
    water_indices = []
    for i in range(len(neighbors)):
        indices = neighbors[i]
        for j in range(len(indices)):
            if indices[j] not in water_indices and indices[j] in traj_solv:
                water_indices.append(indices[j])
    print(len(water_indices))
    return water_indices

def prot_water_dist(traj, prot_res, offset, water_indices):
    from itertools import product
    import mdtraj as md
    import numpy as np

    #Seperate protein residues of interest
    prot_ind = traj_sect(traj, prot_res, offset)

    #Determine all solvent protein residue pairs
    prot_solv = list(product(prot_ind, water_indices))

    #Determine distance b/w all protein residues and all solvent molecules
    [dist, pairs] = md.compute_contacts(traj, contacts=prot_solv, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)

    return dist
