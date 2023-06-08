#!/ usr / bin / env python

#Determine equilibration point for trajectory
#Input: rmsd = BB RMSD for the full trajectory in nm; t_max = length of trajectory in ns) 
def equil_deter(rmsd, t_max, threshold=0.1, per=10, output_ns=True):
    import numpy as np
    
    per = per/100

    #Average rmsd every 200 ps
    int_per_ns = int(len(rmsd)/(t_max))
    rmsd_max = np.zeros(np.round(t_max*5))
    rmsd_min = np.zeros(np.round(t_max*5))
    n=0

    for i in range(len(rmsd_max)):
        k = int(n + (int_per_ns/5))
        rmsd_max[i] = max(rmsd[n:k])
        rmsd_min[i] = min(rmsd[n:k])
        n = k
    #Set count limit
    count_lim = int(len(rmsd)*per)

    #Determine Equilibration time
    time = np.linspace(0, t_max, num = len(rmsd_max))
    eq_time = 'NA'
    count = 0
    for i in range(1, len(rmsd_max)):
        diff = abs(rmsd_max[i] - rmsd_min[i-1])
        if diff < threshold:
            count += 1
        else:
            count = 0
        if count > count_lim:
            eq_time = 5 * (round(time[i-count_lim]/5) + (time[i-count_lim] % 5 > 0)) #round up to nearest 5 ns
            start_i = i-count_lim
            break
    if eq_time == 'NA':
        raise Exception('equilibrium not reched')

    if output_ns == True:
        return eq_time
    else:
        return start_i

#Determine the indices for uncorrelated data
#Input: data = input data
#Output: t_uncorr = indices of the uncorrelated data
def uncorr_ind(data):
    #Import packages
    import ruptures as rpt 
    import numpy as np
    from statistics import stdev
    
    #Raise Error if data trajectory is of length 0
    if len(data) == 0:
        raise Exception('Error! Empty Array Supplied!')

    #Convert data to float
    raw = np.zeros(len(data))
    for i in range(len(data)):
        raw[i] = float(data[i])

    #Apply ruptures to find uncorrelated samples
    model = 'l1'
    algo = rpt.Binseg(model=model, min_size=10, jump=10).fit(raw)
    n = len(raw)
    sigma = stdev(raw)
    #t_uncorr = algo.predict(pen = np.log(n) * sigma**2)
    t_uncorr = algo.predict(epsilon=0.2 * n * sigma ** 2)
    
    #Raise Warning if Ruptures Supplied fewer than 100 uncorrelated samples
    if len(t_uncorr) < 100:
        raise Warning('Low Number of Uncorrelated Samples Found! Review supplied array! n_uncorr = ' + str(len(t_uncorr)))
    
    return t_uncorr

#Sort data to remove correlated samples
#Input: data = full data array, t_uncorr = indices of uncorrelated data
#Output: data_uncorr = data array with correlated samples removed
def uncorr_sort(data, t_uncorr):
    #import packages
    import numpy as np
    
    #Check length of input array against uncorrelated frames list
    if len(data) < t_uncorr[-1]:
        raise Warning('Warning: Supplied array is shorter than expected. Uncorrelated samples not removed!')
        return data
    
    #Check data type and convert if needed or error if not accepted
    if isinstance(data, np.ndarray):
        raw = data
    elif isinstance(data, list):
        #Convert data to float
        raw = np.zeros(len(data))
        for i in range(len(data)):
            raw[i] = float(data[i])
    else:
        raise Exception('Supplied array is not an accepted type (list or np.array)')

    #Reduce to uncorrelated data
    data_uncorr = np.zeros(len(t_uncorr))
    n = 0
    for i in range(len(raw)):
        if i in t_uncorr:
            data_uncorr[n] = raw[i]
            n += 1

    return data_uncorr

#Sort data to remove correlated samples for non-interger data arrays
#Input: data = full data array, t_uncorr = indices of uncorrelated data
#Output: data_uncorr = data array with correlated samples removed
def uncorr_char(data, t_uncorr):
    #Reduce to uncorrelated data
    num=len(t_uncorr)
    data_uncorr = []
    for i in range(len(data)):
        if i in t_uncorr:
            data_uncorr.append(data[i])

    return data_uncorr

#Output: Uncorrelated RMSD values and array of uncorrelated frames
#Input: traj = MDTraj trajectory to compute the RMSD, ref = MDTraj reference structure
#Output:
#rmsd_uncorr = RMSD for each frame at each uncorrelated frame in the trajectory
#t = uncorreated time frames in trajectory
def compute_rmsd(traj, ref, t_uncorr = 'none'):
    import mdtraj as md
    import sys
    #Check that trajectory and reference are of same length
    if traj.n_atoms != ref.n_atoms:
        raise Warning('Reference and trajectory are of different lengths!')
    
    rmsd = md.rmsd(traj, ref, parallel=True, precentered=False)

    if t_uncorr == 'none':
        t = uncorr_ind(rmsd)
    else: 
        t = t_uncorr
    rmsd_uncorr = uncorr_sort(rmsd, t)

    if t_uncorr == 'none':
        return rmsd_uncorr, t
    else:
        return rmsd_uncorr

