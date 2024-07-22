#!/usr/bin/env python

# Importing the necessary libraries for plotting
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

pp = PdfPages('rmsds.pdf')

my_palette = ['#355EE7', '#CB201A', '#1E9965', '#D291BC', '#875C36', '#FEB301']
sns.set_theme(style="ticks")
sns.set_palette(sns.color_palette(my_palette))

# Loading data from the rmsd.xvg file, ignoring lines starting with "@" or "#"
#t1,rmsd_calpha = np.loadtxt("./data/rmsd-calpha.xvg", comments=["@", "#"], unpack=True)
t2,rmsd_prot = np.loadtxt("./data/rmsd-protein.xvg", comments=["@", "#"], unpack=True)
#t3,rmsd_epg = np.loadtxt("./data/rmsd-epg.xvg", comments=["@", "#"], unpack=True)
#t4,rmsd_epg2 = np.loadtxt("./data/rmsd-epg2.xvg", comments=["@", "#"], unpack=True)
#t5,rmsd_epg3 = np.loadtxt("./data/rmsd-epg3.xvg", comments=["@", "#"], unpack=True)
 
# Creating a figure and axis objects for the plot
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)

# Filling the area between the x-values (res) and y-values (rmsf) with a semi-transparent black color
#ax.fill_between(t,rmsd, color="black", linestyle="-", alpha=0.3)

# Plotting the line representing the RMSF values
ax.plot(t2,rmsd_prot, label = "Protein")
#ax.plot(t1,rmsd_calpha, label =r"C$_\alpha$")
#ax.plot(t3,rmsd_epg, label ="EPG1")
#ax.plot(t4,rmsd_epg2, label ="EPG2")
#ax.plot(t5,rmsd_epg3, label ="EPG3")
#plt.legend(ncol = 3)

# Setting labels for the x-axis (time in ps) and y-axis (RMSD value)
ax.set_xlabel("time (ns)")
ax.set_ylabel(r"RMSD (nm)")

# Saving the plot as a PNG image with higher resolution (300 dpi)
#plt.savefig("rmsd.png", format="png", dpi=300)
pp.savefig()
pp.close()

# Displaying the plot
plt.show()
