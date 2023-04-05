#!/ usr / bin / env python
#Import packages
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Equilibration Time from RMSD')
parser.add_argument('-n', required=True, help='File name for BB RMSD (xvg)')
parser.add_argument('-t', required=True, type=int, help='Total time for trajectory (ns)')

#Import Arguments
args = parser.parse_args()
file_name = args.n
t_max = args.t

#Import custom modules
file_path = sys.path[0]
repo_path = file_path.rsplit('/',1)[0]
sys.path.insert(1, repo_path + '/Traj_process')
import process_traj
import load_data

t, rmsd = load_data.col2_float_data('.', file_name, False)

eq_time = process_traj.equil_deter(rmsd, t_max, True)

#Save equilibration time to file
output = open('equilibration_time.txt', 'w')
output.write('Equilibration time = ' + str(eq_time) + ' ns')
