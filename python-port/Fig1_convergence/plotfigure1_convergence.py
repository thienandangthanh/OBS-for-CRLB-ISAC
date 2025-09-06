import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

import sys
import os
# Add parent directory to path to import constants
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from constants import COLORS

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Plot Figure 1 Convergence Analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-file', type=str, default='data_convergence.mat',
                        help='Input MATLAB data file')
    parser.add_argument('--output-file', type=str, default='fig1_convergence_combined.png',
                        help='Output PNG file name')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory')
    parser.add_argument('--save-only', action='store_true',
                        help='Save plot without displaying it')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Resolution for saved plot')
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load the MATLAB data file
    try:
        data = scipy.io.loadmat(args.input_file)
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found.")
        print("Make sure to run propose_SCA.py first to generate the data file.")
        return
    except Exception as e:
        print(f"Error loading data file: {e}")
        return

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

    # Create the plot
    plt.figure(1, figsize=(7, 5))
    plt.xlabel('Number of Iterations')
    plt.ylabel('Objective value')
    plt.grid(True)

    # Plot the three lines
    plt.plot(Z1, color=COLORS[0], linewidth=1.5, linestyle='-', label='δₛ=1, δc=0.05')
    plt.plot(Z2, color=COLORS[1], linewidth=1.5, linestyle='-', label='δₛ=1, δc=0.1')
    plt.plot(Z3, color=COLORS[2], linewidth=1.5, linestyle='-', label='δₛ=1, δc=0.15')

    # Add legend
    plt.legend(loc='lower right')  # 'southeast' in MATLAB corresponds to 'lower right'

    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(args.output_dir, args.output_file)
    plt.savefig(output_path, dpi=args.dpi, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    
    # Display the plot unless save-only mode
    if not args.save_only:
        plt.show()
    else:
        plt.close()

if __name__ == "__main__":
    main()
