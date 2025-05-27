import scipy.io
import matplotlib.pyplot as plt
import numpy as np

import sys
import os
# Add parent directory to path to import constants
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from constants import COLORS

# Load the MATLAB data files
data1 = scipy.io.loadmat('data_WMMSESDR_K.mat')
data2 = scipy.io.loadmat('data_FP_SGDA.mat')
data3 = scipy.io.loadmat('data_SCA_Ns=3M.mat')

# Extract and average SR_all and CRB_all
SR1 = np.mean(data1['SR_all'], axis=0)
SR2 = np.mean(data2['SR_all'], axis=0)
SR3 = np.mean(data3['SR_all'], axis=0)

CRB1 = np.mean(data1['CRB_all'], axis=0)
CRB2 = np.mean(data2['CRB_all'], axis=0)
CRB3 = np.mean(data3['CRB_all'], axis=0)

# Compute the weighted sum
WS1 = 0.25 * SR1 - CRB1
WS2 = 0.25 * SR2 - CRB2
WS3 = 0.25 * SR3 - CRB3

# X-axis labels
x_labels = np.arange(1, 7)  # 1 to 6

# Create figure
plt.figure(figsize=(7, 5))
plt.grid(True)

# Plot SR
line1 = plt.plot(x_labels, SR1, '-.o', linewidth=1.5, color=COLORS[0], label='WMMSE-SDR, SR')[0]
line2 = plt.plot(x_labels, SR2, '-.s', linewidth=1.5, color=COLORS[1], label='FP-SGDA, SR')[0]
line3 = plt.plot(x_labels, SR3, '-.^', linewidth=1.5, color=COLORS[2], label='Algorithm 1, SR')[0]
lines_sr = [line1, line2, line3]

# Plot CRB
line_crb1 = plt.plot(x_labels, CRB1, '--<', linewidth=1.5, color=COLORS[0], label='WMMSE-SDR, CRLB')[0]
line_crb2 = plt.plot(x_labels, CRB2, '-->', linewidth=1.5, color=COLORS[1], label='FP-SGDA, CRLB')[0]
line_crb3 = plt.plot(x_labels, CRB3, '--v', linewidth=1.5, color=COLORS[2], label='Algorithm 1, CRLB')[0]
lines_crb = [line_crb1, line_crb2, line_crb3]

# Plot Weighted Sum (WS)
line_ws1 = plt.plot(x_labels, WS1, '-d', linewidth=1.5, color=COLORS[0], label='WMMSE-SDR, OV')[0]
line_ws2 = plt.plot(x_labels, WS2, '-x', linewidth=1.5, color=COLORS[1], label='FP-SGDA, OV')[0]
line_ws3 = plt.plot(x_labels, WS3, '-p', linewidth=1.5, color=COLORS[2], label='Algorithm 1, OV')[0]
lines_ws = [line_ws1, line_ws2, line_ws3]

# Set labels
plt.xlabel('Number of Communications Users')
plt.ylabel('Metric Values')

# Create three separate legends
# First legend (Sum Rate)
leg1 = plt.legend(lines_sr, [l.get_label() for l in lines_sr], 
                 loc='upper left', title='Sum Rate (nats/Hz)')

# Second legend (CRLB)
ax2 = plt.gca().add_artist(leg1)
leg2 = plt.legend(lines_crb, [l.get_label() for l in lines_crb], 
                 loc='upper right', title='CRLB')

# Third legend (Objective Value)
ax3 = plt.gca().add_artist(leg2)
leg3 = plt.legend(lines_ws, [l.get_label() for l in lines_ws], 
                 loc='center', title='Objective Value')

# Set font properties
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10

# Adjust layout and display
plt.tight_layout()
plt.show()
