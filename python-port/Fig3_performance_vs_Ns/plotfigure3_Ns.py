"""
Figure 3: Performance vs Number of Sensing Streams (Ns)
Python translation of plotfigure3_Ns.m

This script generates plots showing objective function values versus the number
of sensing streams for different algorithms and user configurations.
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
    Load and process data from multiple .mat files.

    Returns:
        dict: Dictionary containing processed data for plotting
    """
    data_dict = {}

    try:
        # Load main data file
        print("Loading main data...")
        main_data = scipy.io.loadmat('data.mat')
        print(f"Main data keys: {list(main_data.keys())}")

        # Extract and process Obj_all data
        # MATLAB: K1=mean(squeeze(Obj_all(2,:,:))); (index 2 in MATLAB = index 1 in Python)
        Obj_all = main_data['Obj_all']
        print(f"Obj_all shape: {Obj_all.shape}")

        # Extract means for K=1, K=2, K=3 (indices 1, 2, 3 in Python)
        data_dict['K1'] = np.mean(np.squeeze(Obj_all[1, :, :]), axis=0)
        data_dict['K2'] = np.mean(np.squeeze(Obj_all[2, :, :]), axis=0)
        data_dict['K3'] = np.mean(np.squeeze(Obj_all[3, :, :]), axis=0)

        # Load sensing data
        print("Loading sensing data...")
        sensing_data = scipy.io.loadmat('data_sensing.mat')
        print(f"Sensing data keys: {list(sensing_data.keys())}")
        data_dict['K0'] = sensing_data['stream_test'].flatten()

        # Load SDR sensing data
        print("Loading SDR sensing data...")
        sdr_sensing_data = scipy.io.loadmat('data_SDR_sensing.mat')
        print(f"SDR sensing data keys: {list(sdr_sensing_data.keys())}")
        data_dict['K0_SDR'] = -sdr_sensing_data['CRB_all'].flatten()

        # Load SDR Ns data and process it
        print("Loading SDR Ns data...")
        sdr_ns_data = scipy.io.loadmat('data_SDR_Ns.mat')
        print(f"SDR Ns data keys: {list(sdr_ns_data.keys())}")
        # Process SDR data similar to MATLAB code
        SR_all = sdr_ns_data['SR_all']
        CRB_all = sdr_ns_data['CRB_all']
        I_out = int(sdr_ns_data.get('I_out', np.array([[100]]))[0, 0])  # Explicitly convert to int
        print(f"SR_all shape: {SR_all.shape}, CRB_all shape: {CRB_all.shape}, I_out: {I_out}")

        # Initialize obj_all_new equivalent to MATLAB: obj_all_new=zeros(3,100,9);
        obj_all_new = np.zeros((3, 100, 9))

        # Process data similar to MATLAB loop
        for k_par in range(1, 301):  # MATLAB: 1:300, Python: range(1, 301)
            channel = ((k_par - 1) % I_out)  # MATLAB: mod((k_par-1),I_out)+1, adjusted for 0-based indexing
            weight = ((k_par - 1) // I_out)   # MATLAB: floor((k_par-1)/I_out)+1, adjusted for 0-based indexing

            if weight < 3 and channel < 100:  # Bounds checking
                obj_all_new[weight, channel, :] = 0.25 * SR_all[:, k_par-1] - CRB_all[:, k_par-1]

        # Extract means for SDR data
        data_dict['K1_SDR'] = np.mean(np.squeeze(obj_all_new[0, :, :]), axis=0)
        data_dict['K2_SDR'] = np.mean(np.squeeze(obj_all_new[1, :, :]), axis=0)
        data_dict['K3_SDR'] = np.full(9, np.mean(np.squeeze(obj_all_new[2, :, 0])))  # Replicate single value 9 times

        # Define indices
        data_dict['index'] = np.arange(0, 9)      # 0:8
        data_dict['index2'] = np.arange(1, 9)     # 1:8

        # Print shapes for debugging
        for key, value in data_dict.items():
            if isinstance(value, np.ndarray):
                print(f"{key} shape: {value.shape}")
            else:
                print(f"{key}: {value}")

    except Exception as e:
        print(f"Error loading data: {e}")
        print("Please ensure all required .mat files are in the current directory:")
        print("- data.mat")
        print("- data_sensing.mat")
        print("- data_SDR_sensing.mat")
        print("- data_SDR_Ns.mat")
        raise

    return data_dict

def create_performance_plot(data_dict):
    """
    Create the performance vs Ns plot.

    Args:
        data_dict (dict): Dictionary containing the processed data
    """
    # Create figure
    plt.figure(figsize=(7, 5))

    # Enable grid
    plt.grid(True)

    # Set axis labels
    plt.xlabel('Number of sensing streams')
    plt.ylabel('Objective value')

    # Plot Algorithm 1 data (solid lines)
    plt.plot(data_dict['index'], data_dict['K1'], 
             color=COLORS[0], linewidth=1.5, linestyle='-', marker='s', 
             markersize=6, label='Algorithm 1, K=1')

    plt.plot(data_dict['index'], data_dict['K2'], 
             color=COLORS[1], linewidth=1.5, linestyle='-', marker='o', 
             markersize=6, label='Algorithm 1, K=2')

    plt.plot(data_dict['index'], data_dict['K3'], 
             color=COLORS[2], linewidth=1.5, linestyle='-', marker='v', 
             markersize=6, label='Algorithm 1, K=3')

    plt.plot(data_dict['index2'], data_dict['K0'], 
             color=COLORS[3], linewidth=1.5, linestyle='-', marker='d', 
             markersize=6, label='Algorithm 1, K=0')

    # Plot WMMSE-SDR data (dashed lines)
    # Note: Using colors from index 0-3 again but with dashed lines to differentiate
    plt.plot(data_dict['index'], data_dict['K1_SDR'], 
             color=COLORS[0], linewidth=1.5, linestyle='--', marker='s', 
             markersize=6, label='WMMSE-SDR, K=1')

    plt.plot(data_dict['index'], data_dict['K2_SDR'], 
             color=COLORS[1], linewidth=1.5, linestyle='--', marker='o', 
             markersize=6, label='WMMSE-SDR, K=2')

    plt.plot(data_dict['index'], data_dict['K3_SDR'], 
             color=COLORS[2], linewidth=1.5, linestyle='--', marker='v', 
             markersize=6, label='WMMSE-SDR, K=3')

    plt.plot(data_dict['index2'], data_dict['K0_SDR'], 
             color=COLORS[3], linewidth=1.5, linestyle='--', marker='d', 
             markersize=6, label='WMMSE-SDR, K=0')

    # Set axis limits
    plt.xlim([0, 8])
    plt.ylim([-12, 0])

    # Add legend
    plt.legend(loc='lower right')

    # Adjust layout and show
    plt.tight_layout()
    plt.show()

def main():
    """
    Main function to execute the plotting script.
    """
    print("Starting Figure 3: Performance vs Number of Sensing Streams")
    print("=" * 60)

    # Load and process data
    data_dict = load_and_process_data()

    # Create the plot
    create_performance_plot(data_dict)

    print("Plot completed successfully!")

if __name__ == "__main__":
    main()
