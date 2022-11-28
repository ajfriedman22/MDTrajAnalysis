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

