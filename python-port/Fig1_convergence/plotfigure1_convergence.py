import scipy.io
import matplotlib.pyplot as plt
import numpy as np

# Load the MATLAB data file
data = scipy.io.loadmat('data_convergence.mat')

# Print the keys to see what variables are available
print("Available variables in the .mat file:")
print(data.keys())

# Access the Con variable (it will be a numpy array of objects)
Con = data['Con']
print(f"Con shape: {Con.shape}")
print(f"Con type: {type(Con)}")

# Extract the third column from each cell
# Con is a 2D array of objects
Con_flat = Con.flatten()
Z1 = Con_flat[0][:, 2] # Python uses 0-based indexing, so column 3 becomes column 2
Z2 = Con_flat[1][:, 2] # Access second element
Z3 = Con_flat[2][:, 2] # Access third element

# Define colors (similar to MATLAB's default color order)
colors = [
    '#0984E3', #Blue
    '#E17055', #Orange
    '#00CEC9', #Cyan
    '#E84393', #Pink
    '#00B894'  #Green
]

# Create the plot
plt.figure(1, figsize=(7, 5))
plt.xlabel('Number of Iterations')
plt.ylabel('Objective value')
plt.grid(True)

# Plot the three lines
plt.plot(Z1, color=colors[0], linewidth=1.5, linestyle='-', label='δₛ=1, δc=0.05')
plt.plot(Z2, color=colors[1], linewidth=1.5, linestyle='-', label='δₛ=1, δc=0.1')
plt.plot(Z3, color=colors[2], linewidth=1.5, linestyle='-', label='δₛ=1, δc=0.15')

# Add legend
plt.legend(loc='lower right')  # 'southeast' in MATLAB corresponds to 'lower right'

# Set font properties
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10

# Adjust layout and display
plt.tight_layout()
plt.show()
