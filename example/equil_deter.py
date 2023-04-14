#!/ usr / bin / env python
#Import packages
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Equilibration Time from RMSD')
parser.add_argument('-n', required=True, help='File name for BB RMSD (xvg)')
parser.add_argument('-t', required=True, type=float, help='Total time for trajectory (ns)')
parser.add_argument('-p', required=True, type=float, help='Percent of trajectory for which RMSD should be stable')
parser.add_argument('-l', required=True, type=float, help='Allowed fluctuations of RMSD at equilibrium(Angstrom)')

#Import Arguments
args = parser.parse_args()
file_name = args.n
t_max = args.t
per = args.p
threshold = args.l

#Import custom modules
<<<<<<< HEAD
file_path = sys.path[0]
repo_path = file_path.rsplit('/',1)[0]
sys.path.insert(1, repo_path + '/Traj_process')
=======
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 
>>>>>>> 10e4a8adc438eba200e5a983c50b9ba2d79715a5
import process_traj

t, rmsd = load_data.col2_float_data('.', file_name, True)

eq_time = process_traj.equil_deter(rmsd, t_max, threshold, per, True)

#Save equilibration time to file
output = open('equilibration_time.txt', 'w')
output.write('Equilibration time = ' + str(eq_time) + ' ns')
=======
eq_time = [int(eq_time)]
np.savetxt('equilibration_time.txt', eq_time)
print('Equilibration analysis completed')
