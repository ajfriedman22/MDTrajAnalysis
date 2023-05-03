#!/ usr / bin / env python
import numpy as np
import argparse
import pandas as pd
#Declare arguments
parser = argparse.ArgumentParser(description = 'Determine the h-bonds which are different b/w two or more conditions')
parser.add_argument('-f', required=True, help='Input file (condition_name dir1 dir2 ...)')

def split_res(res_atom, res_name_ar, res_num_ar, atom_ar):
    import re
    [res, atom] = res_atom.split('-')
    [res_name, res_num, empty] = re.split('(\d+)', res)
    res_name_ar.append(res_name)
    res_num_ar.append(res_num)
    atom_ar.append(atom)

#Import Arguments
args = parser.parse_args()
File_input = args.f

#Open file paths
file_path = open(File_input, 'r').readlines()

#Load h_bonds in each population
hbonds_all, hbonds_all_full = [],[]
hbonds_condition_name = []
for line in file_path:
    line_sep = line.split(' ')
    hbond_line, hbond_line_full = [],[]
    for n in range(1, len(line_sep)):
        path = line_sep[n].strip()
        input_df = pd.read_csv(path + '/Hbond_per.csv')
        for index, row in input_df.iterrows():
            bond_full = row['Donor Residue Name'] + str(row['Donor Residue ID']) + '-' + row['Donor Atom Name'] + ' -- ' + row['Acceptor Residue Name'] + str(row['Acceptor Residue ID']) + '-' + row['Acceptor Atom Name']
            bond = row['Donor Residue Name'] + str(row['Donor Residue ID']) + '-' + row['Acceptor Residue Name'] + str(row['Acceptor Residue ID'])
            if bond not in hbond_line:
                hbond_line.append(bond)
                hbond_line_full.append(bond_full)
    hbonds_all.append(hbond_line)
    hbonds_all_full.append(hbond_line_full)
    hbonds_condition_name.append(line_sep[0])
print('Files Loaded')

#Determine bonds unique to each condition
for i in range(len(hbonds_condition_name)):
    #Set h-bonds both in and not in set condition
    hbonds_condition = hbonds_all[i]
    hbonds_not_condition =[]
    for k in range(len(hbonds_condition_name)):
        if k != i:
            hbonds_not_condition = hbonds_not_condition + hbonds_all[k]

    for j in range(len(hbonds_condition_name)):
        if j != i:
            for n in hbonds_all[j]:
                if n not in hbonds_not_condition and n not in hbonds_condition:
                    hbonds_not_condition.append(n)
    
    #Determine h-bonds only in the set condition
    hbonds_condition_full = hbonds_all_full[i]
    donor_res_name, donor_res_num, donor_atom_name, acceptor_res_name, acceptor_res_num, acceptor_atom_name = [],[],[],[],[],[]
    for n in range(len(hbonds_condition)):
        if hbonds_condition[n] not in hbonds_not_condition:
            [donor, acceptor] = hbonds_condition_full[n].split(' -- ')
            split_res(donor, donor_res_name, donor_res_num, donor_atom_name)
            split_res(acceptor, acceptor_res_name, acceptor_res_num, acceptor_atom_name)
    df = pd.DataFrame({'Donor Residue Name': donor_res_name, 'Donor Residue ID': donor_res_num, 'Donor Atom Name': donor_atom_name, 'Acceptor Residue Name': acceptor_res_name, 'Acceptor Residue ID': acceptor_res_num, 'Acceptor Atom Name': acceptor_atom_name})
    df.to_csv('Hbonds_only_' + str(hbonds_condition_name[i]) + '.csv')

print('Hbond Correlation Analysis Complete')
