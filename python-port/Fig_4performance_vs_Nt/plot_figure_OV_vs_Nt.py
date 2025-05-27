import scipy.io
import matplotlib.pyplot as plt
import numpy as np

import sys
import os
# Add parent directory to path to import constants
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from constants import COLORS

# Load the MATLAB data files
data_wmmsesdr = scipy.io.loadmat('data_WMMSESDR_Nt.mat')
data_fp_sgda = scipy.io.loadmat('data_FP_SGDA.mat')
data_sca_ns3m = scipy.io.loadmat('data_SCA_Ns=3M.mat')
data_sca_ns0 = scipy.io.loadmat('data_SCA_Ns=0.mat')
data_sca_ld_ns3m = scipy.io.loadmat('data_SCA_LD_Ns=3M.mat')

# Extract Time_all from each dataset
T1 = data_wmmsesdr['Time_all']
T2 = data_fp_sgda['Time_all']
T3 = data_sca_ns3m['Time_all']
T4 = data_sca_ns0['Time_all']
T5 = data_sca_ld_ns3m['Time_all']

# Define x-axis values
x = np.array([8, 16, 32])
x2 = np.array([8, 16, 32, 64, 128])

# Calculate mean times
y1 = np.mean(T1, axis=0)  # Assuming Time_all is stored similarly to MATLAB
y2 = np.mean(T2, axis=0)
y3 = np.mean(T3, axis=0)
y4 = np.mean(T4, axis=0)
y5 = np.mean(T5, axis=0)

# Create the plot
plt.figure(figsize=(7, 5))

# Plot lines with log-log scale
lines = []
lines.append(plt.loglog(x, y1, '-s', color=COLORS[0], linewidth=1.5)[0])
lines.append(plt.loglog(x2, y2, '-*', color=COLORS[1], linewidth=1.5)[0])
lines.append(plt.loglog(x2, y3, '-v', color=COLORS[2], linewidth=1.5)[0])
lines.append(plt.loglog(x2, y4, '-^', color=COLORS[3], linewidth=1.5)[0])
lines.append(plt.loglog(x2, y5, '--d', color=COLORS[4], linewidth=1.5)[0])

# Set labels and grid
plt.xlabel('Number of Transmit Antennas')
plt.ylabel('Average CPU Time (sec)')

# Set x-axis limits and ticks
plt.xlim(8, 128)
plt.xticks(x2, x2)

# Enable grid
plt.grid(True)

# Add legend
plt.legend(['WMMSE-SDR', 'FP-SGDA', 'Algorithm 1, Ns=3M', 'Algorithm 1, Ns=0', 'LD Algorithm1, Ns=3M'],
           loc='upper right')

# Set font properties
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10

# Adjust layout and display
plt.tight_layout()
plt.show()
