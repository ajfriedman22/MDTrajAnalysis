def dssp_remove_space(dssp_list, char_replace):
    frame_max,residue = dssp_list.shape #determine the number of frames and residues for which dssp analysis was completed
    dssp_res_mod = []
    for i in range(frame_max): #loop through each residue seperately
        dssp_res = dssp_list[i,:] #seperate all time values for a single residue
        dssp_res_i = [] 
        for j in dssp_res:
            if j == ' ': #in dssp a space designates a loop or irregular element
                dssp_res_i.append(char_replace) #subsitute an l for this designation to prevent issues with reading character values
            else: #if not a space keep the same character designation
                dssp_res_i.append(j)
        dssp_res_mod.append(dssp_res_i)
    return dssp_res_mod

def per_alpha_helix(dssp, time):
    import numpy as np
    from scipy import stats

    char_num = np.arange(2,15,2)
    
    num = 0
    time_tot = int(len(dssp))
    
    per_ind_tot = round(time_tot/25)
    alpha_per = np.zeros([len(char_num), per_ind_tot])

    #Determine helicity for all residues at each time point
    alpha_per_time = np.zeros(time_tot)
    
    #% helicity for each reaidue
    alpha_char = np.zeros(len(char_num))

    for i in dssp:
        alpha = np.zeros(len(char_num))

        char = i
        c = 0
        #Determine DSSP Values
        for n in char_num:
            if char[n]=='H':
                alpha[c] += 1
                alpha_char[c] += 1
            #Every 20 time steps take a running percentage
            if num % 25 == 0 and num != 0 and num < per_ind_tot*25:
                t = int(num/25)
                alpha_per[c][t] = 100 * alpha_char[c] / 25
                alpha_char[c] = 0
            c+=1
        alpha_per_time[num] = 100* sum(alpha)/(len(char_num))
        #Iterate time step
        num += 1
    
    #Determine overall percent for each residue
    alpha_per_mean = np.zeros(len(char_num))
    alpha_per_sem = np.zeros(len(char_num))

    for i in range(len(char_num)):
        alpha_per_mean[i] = np.mean(alpha_per[i][:])
        alpha_per_sem[i] = stats.sem(alpha_per[i][:])

    if time == False:
        return alpha_per, alpha_per_mean, alpha_per_sem
    else:
        return alpha_per_time

def per_helix(dssp):
    import numpy as np
    from scipy import stats

    char_num = np.arange(2,15,2)
    
    num = 0
    time_tot = int(len(dssp))
    
    per_ind_tot = round(time_tot/25)
    struct_per = np.zeros([len(char_num), per_ind_tot])

    #% helicity for each reaidue
    struct_char = np.zeros(len(char_num))

    for i in dssp:
        struct = np.zeros(len(char_num))

        char = i
        c = 0
        #Determine DSSP Values
        for n in char_num:
            if char[n] == 'H' or char[0] == 'G' or char[0] == 'I':
                struct[c] += 1
                struct_char[c] += 1
            #Every 20 time steps take a running percentage
            if num % 25 == 0 and num != 0 and num < per_ind_tot*25:
                t = int(num/25)
                struct_per[c][t] = 100 * struct_char[c] / 25
                struct_char[c] = 0
            c+=1
        #Iterate time step
        num += 1
    #Determine overall percent for each residue
    struct_per_mean = np.zeros(len(char_num))
    struct_per_sem = np.zeros(len(char_num))

    for i in range(len(char_num)):
        struct_per_mean[i] = np.mean(struct_per[i][:])
        struct_per_sem[i] = stats.sem(struct_per[i][:])

    return struct_per, struct_per_mean, struct_per_sem

