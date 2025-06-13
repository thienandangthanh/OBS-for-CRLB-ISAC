#!/usr/bin/env python3
"""
Plotting script for Figure 6 performance vs transmit power analysis.

This script creates publication-quality plots showing the performance 
of the SCA algorithm versus transmit power (Pt).
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import argparse
import os
import sys

# Add parent directory to path for utils
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import db2pow


def plot_performance_vs_Pt(Obj_all, SR_all, CRB_all, I_in=6, 
                           save_plots=False, show_plot=True, output_dir='.'):
    """
    Plot performance metrics versus transmit power.

    Args:
        Obj_all: Objective values array (I_out, I_in)
        SR_all: Sum rate values array (I_out, I_in)  
        CRB_all: CRB trace values array (I_out, I_in)
        I_in: Number of power levels
        save_plots: Whether to save plots to files
        show_plot: Whether to display plots
        output_dir: Directory to save plots
    """

    # Calculate mean performance across channel realizations
    mean_obj = np.mean(Obj_all, axis=0)
    mean_sr = np.mean(SR_all, axis=0)
    mean_crb = np.mean(CRB_all, axis=0)

    # Create transmit power values (in dBm)
    pt_dbm_values = [5 * (i + 1) - 10 for i in range(I_in)]

    # Set up publication-quality plotting parameters
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 12,
        'font.family': 'serif',
        'axes.linewidth': 1.2,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 11,
        'legend.frameon': True,
        'legend.fancybox': True,
        'legend.shadow': True,
        'grid.alpha': 0.3,
        'lines.linewidth': 2,
        'lines.markersize': 8
    })

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Figure 6: SCA Performance vs Transmit Power', fontsize=16, fontweight='bold')

    # Plot 1: Objective Value vs Transmit Power
    ax1 = axes[0, 0]
    ax1.plot(pt_dbm_values, mean_obj, 'b-o', label='Objective Value', linewidth=2, markersize=6)
    ax1.set_xlabel('Transmit Power (dBm)')
    ax1.set_ylabel('Total Objective Value')
    ax1.set_title('(a) Objective Value vs Pt')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Sum Rate vs Transmit Power
    ax2 = axes[0, 1]
    ax2.plot(pt_dbm_values, mean_sr, 'r-s', label='Sum Rate', linewidth=2, markersize=6)
    ax2.set_xlabel('Transmit Power (dBm)')
    ax2.set_ylabel('Sum Rate (nat/s/Hz)')
    ax2.set_title('(b) Sum Rate vs Pt')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Plot 3: CRB Trace vs Transmit Power
    ax3 = axes[1, 0]
    ax3.semilogy(pt_dbm_values, mean_crb, 'g-^', label='CRB Trace', linewidth=2, markersize=6)
    ax3.set_xlabel('Transmit Power (dBm)')
    ax3.set_ylabel('Trace of CRB')
    ax3.set_title('(c) CRB Trace vs Pt')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # Plot 4: Combined normalized performance
    ax4 = axes[1, 1]
    # Normalize values for comparison
    norm_obj = (mean_obj - np.min(mean_obj)) / (np.max(mean_obj) - np.min(mean_obj))
    norm_sr = (mean_sr - np.min(mean_sr)) / (np.max(mean_sr) - np.min(mean_sr))
    norm_crb = (mean_crb - np.min(mean_crb)) / (np.max(mean_crb) - np.min(mean_crb))

    ax4.plot(pt_dbm_values, norm_obj, 'b-o', label='Objective (norm)', linewidth=2, markersize=6)
    ax4.plot(pt_dbm_values, norm_sr, 'r-s', label='Sum Rate (norm)', linewidth=2, markersize=6)
    ax4.plot(pt_dbm_values, 1 - norm_crb, 'g-^', label='1 - CRB (norm)', linewidth=2, markersize=6)
    ax4.set_xlabel('Transmit Power (dBm)')
    ax4.set_ylabel('Normalized Performance')
    ax4.set_title('(d) Normalized Performance Comparison')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    ax4.set_ylim([0, 1.1])

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save plot if requested
    if save_plots:
        plot_filename = os.path.join(output_dir, 'fig6_performance_vs_Pt.png')
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f'Plot saved to {plot_filename}')
    if show_plot:
        plt.show()
    else:
        plt.close()


def plot_detailed_analysis(Obj_all, SR_all, CRB_all, I_in=6, 
                           save_plots=False, show_plot=True, output_dir='.'):
    """
    Create detailed analysis plots including error bars and distribution plots.
    """

    # Calculate statistics
    mean_obj = np.mean(Obj_all, axis=0)
    std_obj = np.std(Obj_all, axis=0)
    mean_sr = np.mean(SR_all, axis=0)
    std_sr = np.std(SR_all, axis=0)
    mean_crb = np.mean(CRB_all, axis=0)
    std_crb = np.std(CRB_all, axis=0)

    # Create transmit power values
    pt_dbm_values = [5 * (i + 1) - 10 for i in range(I_in)]

    # Create detailed figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Figure 6: Detailed Performance Analysis vs Transmit Power', fontsize=16, fontweight='bold')

    # Plot 1: Objective with error bars
    ax1 = axes[0]
    ax1.errorbar(pt_dbm_values, mean_obj, yerr=std_obj, fmt='b-o', 
                 capsize=5, capthick=2, linewidth=2, markersize=6)
    ax1.set_xlabel('Transmit Power (dBm)')
    ax1.set_ylabel('Total Objective Value')
    ax1.set_title('(a) Objective ± Std Dev')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Sum Rate with error bars
    ax2 = axes[1]
    ax2.errorbar(pt_dbm_values, mean_sr, yerr=std_sr, fmt='r-s',
                 capsize=5, capthick=2, linewidth=2, markersize=6)
    ax2.set_xlabel('Transmit Power (dBm)')
    ax2.set_ylabel('Sum Rate (nat/s/Hz)')
    ax2.set_title('(b) Sum Rate ± Std Dev')
    ax2.grid(True, alpha=0.3)

    # Plot 3: CRB with error bars
    ax3 = axes[2]
    ax3.errorbar(pt_dbm_values, mean_crb, yerr=std_crb, fmt='g-^',
                 capsize=5, capthick=2, linewidth=2, markersize=6)
    ax3.set_xlabel('Transmit Power (dBm)')
    ax3.set_ylabel('Trace of CRB')
    ax3.set_title('(c) CRB Trace ± Std Dev')
    ax3.set_yscale('log')
    ax3.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save detailed plot if requested
    if save_plots:
        plot_filename = os.path.join(output_dir, 'fig6_detailed_analysis_vs_Pt.png')
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f'Detailed analysis plot saved to {plot_filename}')

    if show_plot:
        plt.show()
    else:
        plt.close()


def main():
    """Main function for standalone plotting."""

    parser = argparse.ArgumentParser(description='Plot Figure 6 Performance vs Transmit Power')
    parser.add_argument('--data-file', type=str, default='data_SCA_Ns=3M.mat',
                        help='Input data file (.mat format)')
    parser.add_argument('--save-plots', action='store_true',
                        help='Save plots as files instead of displaying')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save output files')
    parser.add_argument('--detailed', action='store_true',
                        help='Create detailed analysis plots with error bars')

    args = parser.parse_args()

    # Load data from .mat file
    try:
        print(f'Loading data from {args.data_file}...')
        data = sio.loadmat(args.data_file)

        Obj_all = data['Obj_all']
        SR_all = data['SR_all']
        CRB_all = data['CRB_all']
        I_in = int(data['I_in'][0, 0]) if 'I_in' in data else 6

        print(f'Data loaded successfully!')
        print(f'Data shapes: Obj_all={Obj_all.shape}, SR_all={SR_all.shape}, CRB_all={CRB_all.shape}')
        print(f'Power levels: {I_in}')

    except FileNotFoundError:
        print(f'Error: Data file {args.data_file} not found!')
        print('Please run the SCA algorithm first to generate data.')
        return
    except Exception as e:
        print(f'Error loading data: {e}')
        return

    # Create main performance plots
    plot_performance_vs_Pt(
        Obj_all, SR_all, CRB_all,
        I_in=I_in,
        save_plots=args.save_plots,
        show_plot=not args.save_plots,
        output_dir=args.output_dir
    )

    # Create detailed analysis if requested
    if args.detailed:
        plot_detailed_analysis(
            Obj_all, SR_all, CRB_all,
            I_in=I_in,
            save_plots=args.save_plots,
            show_plot=not args.save_plots,
            output_dir=args.output_dir
        )

    print('Plotting completed!')


if __name__ == "__main__":
    main()
