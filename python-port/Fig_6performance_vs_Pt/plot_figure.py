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
    parser = argparse.ArgumentParser(description='Plot Figure 6: Performance vs Transmit Power')
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
        data_dict['wmmsesdr'] = scipy.io.loadmat('data_WMMSESDR_pt.mat')
        data_dict['fp_sgda'] = scipy.io.loadmat('data_FP_SGDA.mat')
        data_dict['sca_ns3m'] = scipy.io.loadmat('data_SCA_Ns=3M.mat')

        print("Successfully loaded all data files")

    except Exception as e:
        print(f"Error loading data: {e}")
        print("Please ensure all required .mat files are in the current directory:")
        print("- data_WMMSESDR_pt.mat")
        print("- data_FP_SGDA.mat")
        print("- data_SCA_Ns=3M.mat")
        raise

    return data_dict

def create_performance_plot(data_dict, save_plots=False, output_dir='.'):
    """
    Create the performance vs transmit power plot.

    Args:
        data_dict (dict): Dictionary containing the loaded data
        save_plots (bool, optional): Whether to save the plots. Defaults to False.
        output_dir (str, optional): Directory to save plot. Defaults to current directory.
    """
    # Extract and average SR_all and CRB_all
    SR1 = np.mean(data_dict['wmmsesdr']['SR_all'], axis=0)
    SR2 = np.mean(data_dict['fp_sgda']['SR_all'], axis=0)
    SR3 = np.mean(data_dict['sca_ns3m']['SR_all'], axis=0)

    CRB1 = np.mean(data_dict['wmmsesdr']['CRB_all'], axis=0)
    CRB2 = np.mean(data_dict['fp_sgda']['CRB_all'], axis=0)
    CRB3 = np.mean(data_dict['sca_ns3m']['CRB_all'], axis=0)

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
    plt.xlabel('Transmit Power (dBm)')
    plt.ylabel('Metric Values')

    # Create three separate legends
    # First legend (Sum Rate)
    leg1 = plt.legend(lines_sr, [l.get_label() for l in lines_sr], 
                      loc='upper center', title='Sum Rate')

    # Second legend (CRLB)
    ax2 = plt.gca().add_artist(leg1)
    leg2 = plt.legend(lines_crb, [l.get_label() for l in lines_crb], 
                      loc='upper right', title='CRLB')

    # Third legend (Objective Value)
    ax3 = plt.gca().add_artist(leg2)
    leg3 = plt.legend(lines_ws, [l.get_label() for l in lines_ws], 
                      loc='lower right', title='Objective Value')

    # Adjust layout
    plt.tight_layout()

    # Save plot if save_plots is True
    if save_plots:
        plot_file = os.path.join(output_dir, 'fig6_algorithm_comparison_vs_Pt.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {plot_file}")

    plt.show()

def main():
    """
    Main function to execute the plotting script.
    """
    print("Starting Figure 6: Performance vs Transmit Power")
    print("=" * 70)

    # Parse command line arguments
    args = parse_args()

    # Load data
    data_dict = load_data()

    # Create the plot
    create_performance_plot(data_dict, args.save_plots, args.output_dir)

    print("Plot completed successfully!")

if __name__ == "__main__":
    main()
