
def bond_per(traj_ns, hbonds):
    import mdtraj as md
    import numpy as np

    per = [] #Declare empty array for percentage of time h-bond is formed
    da_distances = md.compute_distances(traj_ns, hbonds[:,[1,2]], periodic=False) #Compute distance between h-bond donor and acceptor
    da_angles = md.compute_angles(traj_ns, hbonds[:,:], periodic=False) #Compute angle between h-bond donor and acceptor
    [num_t, num_h] = np.shape(da_distances) #save values for number of frames(num_t) and number of bonds(num_b) to caculate
    for j in range(num_h): #Loop through all h-bonds
        count = 0 #Initialize count
        for i in range(num_t): #Loop through all frames
            if da_distances[i,j] <= 0.25 and da_angles[i,j] >= 2.094: #If distance between donor and acceptor is less than 2.5A and the angle is greater than 120 degrees or ~ 2.094 radians
                count +=1
        per.append(100*count/num_t) #Percentage of time each h-bond is present in trajectory
    return per

def deter_bond(top, res1, res2, name1, name2):
    import numpy as np
    import mdtraj as md
    import sys

    donor = top.select('resid ' + str(res1) + ' and name ' + str(name1))
    acceptor = top.select('resid ' + str(res2) + ' and name ' + str(name2))
    H = top.select("resid " + str(res1) + " and element H")
    if len(donor) == 0:
        print('Atom ' + str(name1) + ' does not exist in residue ' + str(res1))
        sys.exit()
    elif len(acceptor) == 0:
        print('Atom ' + str(name2) + ' does not exist in residue ' + str(res2))
        sys.exit()
    elif len(H) == 0:
        print('Residue ' + str(res1) + ' contains no H elements')
        sys.exit()
    return donor, acceptor, H

def deter_H(acceptor, H, traj_ns):
    import numpy as np
    import mdtraj as md
    from itertools import product
    import sys

    #Check that multiple acceptor H pairs submitted
    if len(acceptor) == 0 or len(H) == 0:
        print('Error either Acceptors or Hydrogen list is Empty')
        sys.exit()
    #Measure distance between all hydrogens and the acceptor atom
    bond_d = list(product(acceptor, H))
    dist_all = md.compute_distances(traj_ns, bond_d, periodic = False)
    
    #Determine the minimum mean distance
    mean_dist = np.zeros(len(H))
    for j in range(len(H)):
        mean_dist[j] = np.mean(dist_all[:,j])
    #Determine index for minimum distance
    index_min = np.argmin(mean_dist)
        
    #Atom number for hydrogen likely to be involved in bond
    H_min = H[index_min]
    dist = dist_all[:,index_min]

    return H_min, dist

def hbond_read(bond, offset):
    import sys
    import os.path
    
    sys.path.insert(1, '../Traj_process/')
    import data_process 

    line = data_process.split(bond)

    #Determine indicise of dashed
    ind = []
    for j in range(len(line)):
        if line[j] == '-':
            ind.append(j)
    res1 = str(int(data_process.sep_num(data_process.convert(line[0:ind[0]]))) - offset)
    name1 = ''.join([i for i in data_process.convert(line[0:ind[0]]) if not i.isdigit()])
    atom_name1 = data_process.convert(line[ind[0]+1:ind[1]]).strip()
    res2 = str(int(data_process.sep_num(data_process.convert(line[ind[2]+1:ind[3]]))) - offset)
    name2 = ''.join([i for i in data_process.convert(line[ind[2]+1:ind[3]]) if not i.isdigit()])
    atom_name2 = data_process.convert(line[ind[3]+1:]).strip()
    
    return res1, atom_name1, name1, res2, atom_name2, name2

#Seperate residue names and numbers from MDTraj formmetted h-bonds
def process_bond(bond, offset):
    import sys
    sys.path.insert(1, '../Traj_process/')
    import data_process 

    line = data_process.split(bond.strip())
    #Determine indicise of dashed
    ind = []
    for j in range(len(line)):
        if line[j] == '-':
            ind.append(j)
    res1 = str(int(data_process.sep_num(data_process.convert(line[0:ind[0]]))) - offset)
    name1 = data_process.convert(line[ind[0]+1:ind[1]]).strip()
    res2 = str(int(data_process.sep_num(data_process.convert(line[ind[2]+1:ind[3]]))) - offset)
    name2 = data_process.convert(line[ind[3]+1:]).strip()
    return res1, name1, res2, name2
