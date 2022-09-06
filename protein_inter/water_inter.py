
def prot_water_dist(traj, prot_res, offset):
    from itertools import product

    #Seperate protein residues of interest
    res1 = prot_res[0] - offset
    res2 = prot_res[1] - offset
    top = traj.topology()
    traj_sect = top.select(str(res1) + ' <= resid and resid <= ' + str(res2)) #Select only atoms in the given section
    
    #Seperate solvent molecules
    traj_solv = top.select('resname SOL')
    
    #Determine all solvent protein residue pairs
    prot_solv = list(product(traj_select, traj_solv))

    #Determine distance b/w all protein residues and all solvent molecules
    [dist, pairs] = md.compute_contacts(traj, contacts=prot_solv, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)

    return dist
