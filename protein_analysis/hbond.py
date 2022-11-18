
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

def deter_bond(top, res1, res2, name1, name2, i):
    import numpy as np
    import mdtraj as md

    bond = np.zeros(3)
    donor = top.select('resid ' + str(res1[i]) + ' and name ' + str(name1[i]))
    acceptor = top.select('resid ' + str(res2[i]) + ' and name ' + str(name2[i]))
    H = top.select("resid " + str(res1[i]) + " and element H")
    return donor, acceptor, H

def deter_H(acceptor, H, traj_ns):
    import numpy as np
    import mdtraj as md
    from itertools import product

    #easure distance between all hydrogens and the acceptor atom
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

