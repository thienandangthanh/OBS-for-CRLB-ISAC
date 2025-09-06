"""
Figure 2: Trade-off Region Plot
Python translation of plotfigure2_tradeoff_region.m

This script generates the trade-off region plot showing the relationship between
CRLB (Cramér-Rao Lower Bound) and Sum Rate for different algorithms.
"""

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Import only COLORS constant
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from constants import COLORS

def load_and_process_data():
    """
    Load data from the three .mat files and compute mean values.

    Returns:
        dict: Dictionary containing CRB and SR data for each algorithm
    """
    data_dict = {}

    try:
        # Load WMMSE-SDR data
        print("Loading WMMSE-SDR data...")
        wmmse_data = scipy.io.loadmat('data_WMMSESDR_to4.mat')
        print(f"WMMSE data keys: {list(wmmse_data.keys())}")

        # Extract CRB_all and SR_all, compute means
        data_dict['wmmse'] = {
            'CRB': np.mean(wmmse_data['CRB_all'], axis=0),
            'SR': np.mean(wmmse_data['SR_all'], axis=0)
        }

        # Load Proposed SCA data
        print("Loading Proposed SCA data...")
        sca_data = scipy.io.loadmat('data_proposed_SCA.mat')
        print(f"SCA data keys: {list(sca_data.keys())}")

        data_dict['sca'] = {
            'CRB': np.mean(sca_data['CRB_all'], axis=0),
            'SR': np.mean(sca_data['SR_all'], axis=0)
        }

        # Load FP-SGDA data
        print("Loading FP-SGDA data...")
        fp_data = scipy.io.loadmat('data_FP_SGDA.mat')
        print(f"FP data keys: {list(fp_data.keys())}")

        data_dict['fp'] = {
            'CRB': np.mean(fp_data['CRB_all'], axis=0),
            'SR': np.mean(fp_data['SR_all'], axis=0)
        }

        # Print shapes for debugging
        for key, value in data_dict.items():
            print(f"{key} - CRB shape: {value['CRB'].shape}, SR shape: {value['SR'].shape}")

    except Exception as e:
        print(f"Error loading data: {e}")
        print("Please ensure all required .mat files are in the current directory:")
        print("- data_WMMSESDR_to4.mat")
        print("- data_proposed_SCA.mat") 
        print("- data_FP_SGDA.mat")
        raise

    return data_dict

def create_tradeoff_plot(data_dict):
    """
    Create the trade-off region plot.

    Args:
        data_dict (dict): Dictionary containing the processed data
    """
    # Create figure
    plt.figure(figsize=(7, 5))

    # Enable grid
    plt.grid(True)

    # Set axis labels
    plt.xlabel('Trace of the Inverse of the FIM')
    plt.ylabel('Sum Rate (nats/Hz)')

    # Plot data for each algorithm
    algorithms = ['wmmse', 'fp', 'sca']
    markers = ['v', 'd', 's']  # Triangle down, diamond, square
    linestyles = ['-', '-', '--']  # Solid, solid, dashed
    labels = ['WMMSE-SDR', 'FP-SGDA', 'Algorithm 1']

    for i, (alg, marker, linestyle, label) in enumerate(zip(algorithms, markers, linestyles, labels)):
        crb_data = data_dict[alg]['CRB']
        sr_data = data_dict[alg]['SR']

        # Create marker indices (every 2nd point up to 60, similar to MATLAB [1:2:60])
        marker_indices = list(range(0, min(60, len(crb_data)), 2))

        # Plot the line using COLORS constant
        plt.plot(crb_data, sr_data, 
                 color=COLORS[i],
                 linewidth=1.5,
                 linestyle=linestyle,
                 marker=marker,
                 markevery=marker_indices,
                 markersize=6,
                 label=label)

    # Add legend
    plt.legend(loc='lower right')

    # Adjust layout and show
    plt.tight_layout()
    plt.show()

def main():
    """
    Main function to execute the plotting script.
    """
    print("Starting Figure 2: Trade-off Region Plot")
    print("=" * 50)

    # Load and process data
    data_dict = load_and_process_data()

    # Create the plot
    create_tradeoff_plot(data_dict)

    print("Plot completed successfully!")

if __name__ == "__main__":
    main()
