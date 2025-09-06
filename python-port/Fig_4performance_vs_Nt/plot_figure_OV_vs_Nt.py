import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

import sys
# Add parent directory to path to import constants
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from constants import COLORS

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Plot Figure 4: Average CPU Time vs Number of Transmit Antennas')
    parser.add_argument('--save-plots', action='store_true', help='Save plot as files')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save output files')
    return parser.parse_args()

def load_data():
    """
    Load data from MATLAB .mat files.

    Returns:
        dict: Dictionary containing loaded data
    """
    data_dict = {}

    try:
        # Load the MATLAB data files
        data_dict['wmmsesdr'] = scipy.io.loadmat('data_WMMSESDR_Nt.mat')
        data_dict['fp_sgda'] = scipy.io.loadmat('data_FP_SGDA.mat')
        data_dict['sca_ns3m'] = scipy.io.loadmat('data_SCA_Ns=3M.mat')
        data_dict['sca_ns0'] = scipy.io.loadmat('data_SCA_Ns=0.mat')
        data_dict['sca_ld_ns3m'] = scipy.io.loadmat('data_SCA_LD_Ns=3M.mat')

        print("Successfully loaded all data files")

    except Exception as e:
        print(f"Error loading data: {e}")
        print("Please ensure all required .mat files are in the current directory:")
        print("- data_WMMSESDR_Nt.mat")
        print("- data_FP_SGDA.mat")
        print("- data_SCA_Ns=3M.mat")
        print("- data_SCA_Ns=0.mat")
        print("- data_SCA_LD_Ns=3M.mat")
        raise

    return data_dict

def create_cpu_time_plot(data_dict, save_plots=False, output_dir='.'):
    """
    Create the CPU time vs number of antennas plot.

    Args:
        data_dict (dict): Dictionary containing the loaded data
        save_plots (bool, optional): Whether to save the plots. Defaults to False.
        output_dir (str, optional): Directory to save plot. Defaults to current directory.
    """
    # Extract Time_all from each dataset
    T1 = data_dict['wmmsesdr']['Time_all']
    T2 = data_dict['fp_sgda']['Time_all']
    T3 = data_dict['sca_ns3m']['Time_all']
    T4 = data_dict['sca_ns0']['Time_all']
    T5 = data_dict['sca_ld_ns3m']['Time_all']

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

    # Adjust layout
    plt.tight_layout()

    # Save plot if save_plots is True
    if save_plots:
        plot_file = os.path.join(output_dir, 'plot_figure_OV_vs_Nt.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {plot_file}")

    plt.show()

def main():
    """
    Main function to execute the plotting script.
    """
    print("Starting Figure 4: Average CPU Time vs Number of Transmit Antennas")
    print("=" * 70)

    # Parse command line arguments
    args = parse_args()

    # Load data
    data_dict = load_data()

    # Create the plot
    create_cpu_time_plot(data_dict, args.save_plots, args.output_dir)

    print("Plot completed successfully!")

if __name__ == "__main__":
    main()
