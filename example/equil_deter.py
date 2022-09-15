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
sys.path.insert(1,'/home/anika/Documents/code/MDTrajAnalysis/Traj_process')
import process_traj
import load_data

t, rmsd = load_data.col2_float_data('.', file_name, False)

eq_time = process_traj.equil_deter(rmsd, t_max, True)

print(eq_time)

