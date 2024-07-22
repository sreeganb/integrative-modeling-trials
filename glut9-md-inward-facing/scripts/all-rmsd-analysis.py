#!/usr/bin/env python

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import pandas as pd
import matplotlib.pyplot as plt

# Load the .gro and .xtc files
u = mda.Universe("step6.8_equilibration.gro", "step7_production_rep2_run1.xtc")

# Select the C-alpha atoms and the whole protein
calphas = u.select_atoms("name CA")
protein = u.select_atoms("protein")
urc = u.select_atoms("resname URC")

# Create RMSD objects
R_calpha = rms.RMSD(calphas, calphas, select='name CA', groupselections=["protein and name CA"])
R_protein = rms.RMSD(protein, protein, select='protein')
R_urc = rms.RMSD(urc, urc, select='resname URC')

# Run the RMSD calculations
R_calpha.run()
R_protein.run()
R_urc.run()

# Get the RMSD results
rmsd_results_calpha = R_calpha.rmsd
rmsd_results_protein = R_protein.rmsd
rmsd_results_urc = R_urc.rmsd

# Create DataFrames for the RMSD results
df_calpha = pd.DataFrame(rmsd_results_calpha, columns=["Frame", "Time (ps)", "RMSD", "rmsd2"])
df_protein = pd.DataFrame(rmsd_results_protein, columns=["Frame", "Time (ps)", "RMSD"])
df_urc = pd.DataFrame(rmsd_results_urc, columns=["Frame", "Time (ps)", "RMSD"])

# Save the RMSD results to CSV files
df_calpha.to_csv("rmsd_calpha_results.csv", index=False)
df_protein.to_csv("rmsd_protein_results.csv", index=False)
df_urc.to_csv("rmsd_urc_results.csv", index=False)

# Plot the RMSD data
plt.figure(figsize=(12, 6))

# Plot for C-alpha RMSD
plt.subplot(1, 3, 1)
plt.plot(df_calpha["Time (ps)"], df_calpha["RMSD"], marker='o', linestyle='-', color='b')
plt.xlabel("Time (ps)")
plt.ylabel("RMSD")
plt.title("RMSD of C-alpha Atoms Over Time")
plt.grid(True)

# Plot for whole protein RMSD
plt.subplot(1, 3, 2)
plt.plot(df_protein["Time (ps)"], df_protein["RMSD"], marker='o', linestyle='-', color='r')
plt.xlabel("Time (ps)")
plt.ylabel("RMSD")
plt.title("RMSD of Whole Protein Over Time")
plt.grid(True)

#plot for URC RMSD
plt.subplot(1, 3, 3)
plt.plot(df_urc["Time (ps)"], df_urc["RMSD"], marker='o', linestyle='-', color='g')
plt.xlabel("Time (ps)")
plt.ylabel("RMSD")
plt.title("RMSD of URC Over Time")
plt.grid(True)

plt.tight_layout()
plt.savefig("rmsd_rep2.pdf")
plt.show()

print("RMSD results saved to 'rmsd_calpha_results.csv' and 'rmsd_protein_results.csv'. Plot saved to 'rmsd_plots.png'.")
