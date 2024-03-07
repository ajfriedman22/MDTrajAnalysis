import mdtraj as md
import subprocess
import pandas as pd
import argparse
import pymbar
import numpy as np
import matplotlib.pyplot as plt

#Declare arguments
parser = argparse.ArgumentParser(description = 'Process trajectory and return equilibration time, correlation time, and equilibrated trajecotry')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g

#Load trajectory
traj = md.load(File_traj, File_gro)
traj = traj.remove_solvent()

#Compute RMSD on heavy atoms
rmsd = md.rmsd(traj, atom_indices=traj.topology.select('element != H'))

#Determine equilibration point
[t, g, Neff_max] = pymbar.timeseries.detect_equilibration(rmsd, nskip=5)
print(t)
rmsd_equil = rmsd[t:]

#Determine auto-correlation time from equilibrated data
df = pd.series(rmsd_equil)
auto_corr = np.zeros(len(rmsd_equil)/2-1)
act = None
for n in range(1, len(rmsd_equil)/2):
    auto_corr[n-1] = df.autocorr(lag=n)
    if auto_corr[n-2] > 0 and auto_corr[n-1] < 0:
        act = n

plt.figure()
plt.scatter(auto_corr)
plt.xlabel('lag')
plt.ylabel('Auto Correlation')
plt.savefig('autocorr.png')
print(act)

#grompp = subprocess.run(
#            [
#                self.gmx_executable,
#                "grompp",
#                "-f",
#                mdp_file,
#                "-p",
#                self.parametrized_system,
#                "-c",
#                self.initial_coords_file,
#                "-maxwarn",
#                "2",
#                "-o",
#                str(self.output_prefix),
#            ]
#        )