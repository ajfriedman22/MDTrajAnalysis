#!/ usr / bin / env python

#Determine equilibration point for trajectory
#Input: rmsd = BB RMSD for the full trajectory in nm; t_max = length of trajectory in ns) 
def equil_deter(rmsd, t_max, output_ns):
    import numpy as np

    #Average rmsd every 200 ps
    int_per_ns = int(len(rmsd)/(t_max))
    rmsd_max = np.zeros(t_max*5)
    rmsd_min = np.zeros(t_max*5)
    n=0

    for i in range(len(rmsd_max)):
        k = int(n + (int_per_ns/5))
        rmsd_max[i] = max(rmsd[n:k])
        rmsd_min[i] = min(rmsd[n:k])
        n = k

    #Determine Equilibration time
    time = np.linspace(0, t_max, num = len(rmsd_max))
    eq_time = 'NA'
    count = 0
    for i in range(1, len(rmsd_max)):
        diff = abs(rmsd_max[i] - rmsd_min[i-1])
        if diff < 0.1:
            count += 1
        else:
            count = 0
        if count > 50:
            eq_time = 5 * (round(time[i-50]/5) + (time[i-50] % 5 > 0)) #round up to nearest 5 ns
            start_i = i-50
            break
    if eq_time == 'NA':
        print('equilibrium not reched')

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

    #Convert data to float
    raw = np.zeros(len(data))
    for i in range(len(data)):
        raw[i] = float(data[i])

    #Apply ruptures to find uncorrelated samples
    model = 'l1'
    algo = rpt.Binseg(model=model).fit(raw)
    n = len(raw)
    sigma = stdev(raw)
    t_uncorr = algo.predict(pen = np.log(n) * sigma**2)

    return t_uncorr

#Sort data to remove correlated samples
#Input: data = full data array, t_uncorr = indices of uncorrelated data
#Output: data_uncorr = data array with correlated samples removed
def uncorr_sort(data, t_uncorr):
    #import packages
    import numpy as np

    #Convert data to float
    raw = np.zeros(len(data))
    for i in range(len(data)):
        raw[i] = float(data[i])

    #Reduce to uncorrelated data
    num=len(t_uncorr)
    data_uncorr = np.zeros(num)
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
def compute_rmsd(traj, ref):
    import mdtraj as md
    import sys

    rmsd = md.rmsd(traj, ref, parallel=True, precentered=False)
    t = uncorr_ind(rmsd)
    rmsd_uncorr = uncorr_sort(rmsd, t)
    return rmsd_uncorr, t

