#!/ usr / bin / env python
import numpy as np
import argparse
import pandas as pd
#Declare arguments
parser = argparse.ArgumentParser(description = 'Determine the h-bonds which are different b/w two or more conditions')
parser.add_argument('-f', required=True, help='Input file (condition_name dir1 dir2 ...)')

#Import Arguments
args = parser.parse_args()
File_input = args.f

#Open file paths
file_path = open(File_input, 'r').readlines()

#Load h_bonds in each population
hbonds_all = []
hbonds_condition_name = []
for line in file_path:
    line_sep = line.split(' ')
    hbond_line = []
    for n in range(1, len(line_sep)):
        path = line_sep[n].strip()
        input_df = pd.read_csv(path + '/Hbonds_per.csv')
        for index, row in input_df.iterrows():
            bond = row['Donor Residue Name'] + row['Donor Residue Num'] + row['Donor Atom Num'] + '-' + row['Acceptor Reisude Name'] + row['Acceptor Residue Num'] + row['Acceptor Atom Num']]
            if bond not in hbond_line:
                hbond_line.append(bond)
    hbonds_all.append(hbond_line)
    hbonds_condition_name.append(line_sep[0])
print('Files Loaded')

#Determine bonds unique to each condition
for i in range(len(hbonds_condition_name)):
    #Set h-bonds both in and not in set condition
    hbonds_condition = hbonds_all[i]
    hbonds_not_condition = []
    for j in range(len(hbonds_condition_name)):
        if j != i:
            for n in hbonds_all[j]:
                if n not in hbonds_not_condition and n not in hbonds_condition:
                    hbonds_not_condition.append(n)
    
    #Determine h-bonds only in the set condition
    output_hbonds = open('Hbonds_only_' + str(hbonds_condition_name[i]) + '.txt', 'w')
    for n in hbonds_condition:
        if n not in hbonds_not_condition:
            output_hbonds.write(n)
print('Hbond Correlation Analysis Complete')
