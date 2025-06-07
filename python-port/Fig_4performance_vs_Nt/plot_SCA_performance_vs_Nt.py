#!/usr/bin/env python3
"""
Plotting module for Figure 4 SCA performance vs number of antennas (Nt).

This module provides functions to visualize the results from the proposed_SCA.py
implementation, including performance metrics vs antenna count and convergence plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import sys
from typing import Optional, List, Tuple

# Add parent directory to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

try:
    from constants import COLORS
except ImportError:
    # Fallback colors if constants module not available
    COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


def plot_performance_vs_Nt(
    SR_all: np.ndarray,
    CRB_all: np.ndarray,
    Time_all: np.ndarray,
    Obj_all: np.ndarray,
    antenna_sizes: Optional[List[int]] = None,
    save_plots: bool = False,
    show_plot: bool = True,
    output_dir: str = '.',
    figure_title: str = 'Performance vs Number of Transmit Antennas'
) -> None:
    """
    Plot performance metrics vs number of transmit antennas.
    
    Args:
        SR_all: Sum rate results array (I_out, I_in)
        CRB_all: CRB trace results array (I_out, I_in)
        Time_all: Computation time results array (I_out, I_in)
        Obj_all: Objective function results array (I_out, I_in)
        antenna_sizes: List of antenna counts. If None, uses [8,16,32,64,128]
        save_plots: Whether to save plots to files
        show_plot: Whether to display plots
        output_dir: Directory to save plots
        figure_title: Title for the main figure
    """
    
    if antenna_sizes is None:
        antenna_sizes = [8, 16, 32, 64, 128]
    
    # Calculate mean and std for each antenna configuration
    mean_sr = np.mean(SR_all, axis=0)
    std_sr = np.std(SR_all, axis=0)
    
    mean_crb = np.mean(CRB_all, axis=0)
    std_crb = np.std(CRB_all, axis=0)
    
    mean_time = np.mean(Time_all, axis=0)
    std_time = np.std(Time_all, axis=0)
    
    mean_obj = np.mean(Obj_all, axis=0)
    std_obj = np.std(Obj_all, axis=0)
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Objective Function vs Nt
    ax1.errorbar(antenna_sizes, mean_obj, yerr=std_obj, 
                 marker='o', linewidth=2, capsize=5, color=COLORS[0])
    ax1.set_xlabel('Number of Transmit Antennas')
    ax1.set_ylabel('Objective Value')
    ax1.set_title('Objective Function vs Nt')
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log', base=2)
    ax1.set_xticks(antenna_sizes)
    ax1.set_xticklabels(antenna_sizes)
    
    # Plot 2: Sum Rate vs Nt
    ax2.errorbar(antenna_sizes, mean_sr, yerr=std_sr,
                 marker='s', linewidth=2, capsize=5, color=COLORS[1])
    ax2.set_xlabel('Number of Transmit Antennas')
    ax2.set_ylabel('Sum Rate (nat/s/Hz)')
    ax2.set_title('Sum Rate vs Nt')
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log', base=2)
    ax2.set_xticks(antenna_sizes)
    ax2.set_xticklabels(antenna_sizes)
    
    # Plot 3: CRB Trace vs Nt (log scale)
    ax3.errorbar(antenna_sizes, mean_crb, yerr=std_crb,
                 marker='^', linewidth=2, capsize=5, color=COLORS[2])
    ax3.set_xlabel('Number of Transmit Antennas')
    ax3.set_ylabel('CRB Trace')
    ax3.set_title('CRB Trace vs Nt')
    ax3.set_yscale('log')
    ax3.grid(True, alpha=0.3)
    ax3.set_xscale('log', base=2)
    ax3.set_xticks(antenna_sizes)
    ax3.set_xticklabels(antenna_sizes)
    
    # Plot 4: Computation Time vs Nt (log-log scale)
    ax4.errorbar(antenna_sizes, mean_time, yerr=std_time,
                 marker='d', linewidth=2, capsize=5, color=COLORS[3])
    ax4.set_xlabel('Number of Transmit Antennas')
    ax4.set_ylabel('Average CPU Time (sec)')
    ax4.set_title('Computation Time vs Nt')
    ax4.set_yscale('log')
    ax4.set_xscale('log', base=2)
    ax4.grid(True, alpha=0.3)
    ax4.set_xticks(antenna_sizes)
    ax4.set_xticklabels(antenna_sizes)
    
    plt.suptitle(figure_title, fontsize=16)
    plt.tight_layout()
    
    if save_plots:
        plot_file = os.path.join(output_dir, 'fig4_performance_vs_Nt.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Performance plot saved to: {plot_file}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()


def plot_convergence_analysis(
    Con: np.ndarray,
    save_plots: bool = False,
    show_plot: bool = True,
    output_dir: str = '.',
    scenario_info: str = 'Last Scenario'
) -> None:
    """
    Plot convergence analysis for SCA algorithm.
    
    Args:
        Con: Convergence data array (iterations, 3) with columns [sum_rate, -trace_inv_fim, obj]
        save_plots: Whether to save plots to files
        show_plot: Whether to display plots
        output_dir: Directory to save plots
        scenario_info: Information about the scenario being plotted
    """
    
    # Create figure with subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    
    iterations = range(1, len(Con) + 1)
    
    # Plot 1: Sum Rate Convergence
    ax1.plot(iterations, Con[:, 0], color=COLORS[0], linewidth=2)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Sum Rate (nat/s/Hz)')
    ax1.set_title(f'Sum Rate Convergence ({scenario_info})')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Sensing Performance (Trace of inverse FIM)
    ax2.plot(iterations, Con[:, 1], color=COLORS[1], linewidth=2)
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Trace of inverse FIM')
    ax2.set_title(f'Sensing Performance Convergence ({scenario_info})')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Objective Function
    ax3.plot(iterations, Con[:, 2], color=COLORS[2], linewidth=2)
    ax3.set_xlabel('Iteration')
    ax3.set_ylabel('Total Objective Value')
    ax3.set_title(f'Objective Function Convergence ({scenario_info})')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plots:
        plot_file = os.path.join(output_dir, 'fig4_convergence_analysis.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Convergence plot saved to: {plot_file}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()


def plot_comparison_with_other_algorithms(
    our_time_data: np.ndarray,
    antenna_sizes: Optional[List[int]] = None,
    data_dir: str = '.',
    save_plots: bool = False,
    show_plot: bool = True,
    output_dir: str = '.'
) -> None:
    """
    Plot comparison with other algorithms using existing data files.
    
    Args:
        our_time_data: Our algorithm's computation time data (I_out, I_in)
        antenna_sizes: List of antenna counts
        data_dir: Directory containing the MATLAB data files
        save_plots: Whether to save plots to files
        show_plot: Whether to display plots
        output_dir: Directory to save plots
    """
    
    if antenna_sizes is None:
        antenna_sizes = [8, 16, 32, 64, 128]
    
    # Try to load existing algorithm data
    algorithm_data = {}
    
    data_files = {
        'WMMSE-SDR': 'data_WMMSESDR_Nt.mat',
        'FP-SGDA': 'data_FP_SGDA.mat',
        'SCA (Ns=3M)': 'data_SCA_Ns=3M.mat',
        'SCA (Ns=0)': 'data_SCA_Ns=0.mat',
        'LD SCA (Ns=3M)': 'data_SCA_LD_Ns=3M.mat'
    }
    
    for name, filename in data_files.items():
        try:
            filepath = os.path.join(data_dir, filename)
            if os.path.exists(filepath):
                data = sio.loadmat(filepath)
                if 'Time_all' in data:
                    algorithm_data[name] = np.mean(data['Time_all'], axis=0)
                    print(f"Loaded data for {name}")
        except Exception as e:
            print(f"Could not load {name} data: {e}")
    
    # Add our algorithm data
    our_mean_time = np.mean(our_time_data, axis=0)
    algorithm_data['Our SCA (Fig4)'] = our_mean_time
    
    if not algorithm_data:
        print("No algorithm data available for comparison")
        return
    
    # Create comparison plot
    plt.figure(figsize=(10, 7))
    
    # Different antenna sizes for different algorithms
    x_wmmse = [8, 16, 32]  # WMMSE-SDR typically has fewer points
    x_full = antenna_sizes[:len(antenna_sizes)]  # Full range
    
    for i, (name, times) in enumerate(algorithm_data.items()):
        if name == 'WMMSE-SDR' and len(times) == 3:
            x_data = x_wmmse
        else:
            x_data = x_full[:len(times)]
        
        if len(times) > 0:
            if 'WMMSE' in name:
                plt.loglog(x_data, times, '-s', color=COLORS[i % len(COLORS)], 
                          linewidth=1.5, label=name)
            elif 'FP-SGDA' in name:
                plt.loglog(x_data, times, '-*', color=COLORS[i % len(COLORS)], 
                          linewidth=1.5, label=name)
            elif 'SCA (Ns=3M)' in name:
                plt.loglog(x_data, times, '-v', color=COLORS[i % len(COLORS)], 
                          linewidth=1.5, label=name)
            elif 'SCA (Ns=0)' in name:
                plt.loglog(x_data, times, '-^', color=COLORS[i % len(COLORS)], 
                          linewidth=1.5, label=name)
            elif 'LD' in name:
                plt.loglog(x_data, times, '--d', color=COLORS[i % len(COLORS)], 
                          linewidth=1.5, label=name)
            else:
                plt.loglog(x_data, times, '-o', color=COLORS[i % len(COLORS)], 
                          linewidth=1.5, label=name)
    
    plt.xlabel('Number of Transmit Antennas')
    plt.ylabel('Average CPU Time (sec)')
    plt.title('Computation Time Comparison vs Number of Antennas')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper left')
    plt.xlim(8, 128)
    
    # Set x-axis ticks
    plt.xticks(antenna_sizes, antenna_sizes)
    
    plt.tight_layout()
    
    if save_plots:
        plot_file = os.path.join(output_dir, 'fig4_algorithm_comparison.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Comparison plot saved to: {plot_file}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()


def load_and_plot_from_mat_file(
    mat_file: str,
    save_plots: bool = False,
    show_plot: bool = True,
    output_dir: str = '.'
) -> None:
    """
    Load results from .mat file and create all plots.
    
    Args:
        mat_file: Path to the .mat file containing results
        save_plots: Whether to save plots to files
        show_plot: Whether to display plots
        output_dir: Directory to save plots
    """
    
    try:
        data = sio.loadmat(mat_file)
        
        # Extract data arrays
        SR_all = data['SR_all']
        CRB_all = data['CRB_all']
        Time_all = data['Time_all']
        Obj_all = data['Obj_all']
        
        print(f"Loaded data from {mat_file}")
        print(f"Data shape: {SR_all.shape}")
        
        # Create all plots
        plot_performance_vs_Nt(
            SR_all, CRB_all, Time_all, Obj_all,
            save_plots=save_plots, show_plot=show_plot, output_dir=output_dir
        )
        
        # If data directory exists, create comparison plot
        data_dir = os.path.dirname(mat_file)
        plot_comparison_with_other_algorithms(
            Time_all, data_dir=data_dir,
            save_plots=save_plots, show_plot=show_plot, output_dir=output_dir
        )
        
    except Exception as e:
        print(f"Error loading data from {mat_file}: {e}")


def main():
    """Main function for standalone plotting."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Plot Figure 4 SCA performance vs Nt results')
    parser.add_argument('--mat-file', type=str, 
                        default='data_proposed_SCA_fig4.mat',
                        help='Path to .mat file with results')
    parser.add_argument('--save-plots', action='store_true',
                        help='Save plots as files')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save plots')
    
    args = parser.parse_args()
    
    if os.path.exists(args.mat_file):
        load_and_plot_from_mat_file(
            args.mat_file,
            save_plots=args.save_plots,
            show_plot=not args.save_plots,
            output_dir=args.output_dir
        )
    else:
        print(f"Data file {args.mat_file} not found.")
        print("Run proposed_SCA.py first to generate results.")


if __name__ == "__main__":
    main() 