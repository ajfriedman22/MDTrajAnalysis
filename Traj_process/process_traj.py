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
        if count > 100:
            eq_time = 5 * (round(time[i-50]/5) + (time[i-50] % 5 > 0)) #round up to nearest 5 ns
            start_i = i-50
            break
    if eq_time == 'NA':
        print('equilibrium not reched')

    if output_ns == True:
        return eq_time
    else:
        return start_i
