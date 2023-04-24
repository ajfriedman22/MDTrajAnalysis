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
parser.add_argument('-t', required=True, type=int, help='Total time for trajectory (ns)')
parser.add_argument('-p', required=True, type=float, help='Percent of trajectory for which RMSD should be stable')
parser.add_argument('-l', required=True, type=float, help='Allowed fluctuations of RMSD at equilibrium(Angstrom)')

#Import Arguments
args = parser.parse_args()
file_name = args.n
t_max = args.t
per = args.p
threshold = args.l

#Import custom modules
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 
import process_traj

t, rmsd = load_data.col2_float_data('.', file_name, True)

eq_time = process_traj.equil_deter(rmsd, t_max, threshold, per, True)

#Save equilibration time to file
eq_time = [int(eq_time)]
np.savetxt('equilibration_time.txt', eq_time)
print('Equilibration analysis completed')
