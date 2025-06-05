#!/usr/bin/env python3
"""
Plot performance vs number of sensing streams (Ns) for the OBS-for-CRLB-ISAC system.
This module provides both a function that can be imported and a script interface.
"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import os
import sys
import argparse
# Add parent directory to path to import constants
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from constants import COLORS


def plot_performance_vs_Ns(Obj_all, Ns_all=None, save_plots=False, show_plot=True, output_dir='.'):
    """
    Plot performance vs number of sensing streams (Ns).

    Args:
        Obj_all (np.ndarray): Array of objective values with shape (4, I_out, I_in)
        Ns_all (list/array, optional): Array of Ns values. Defaults to [0,1,2,3,4,5,6,7,8]
        save_plots (bool, optional): Whether to save the plots. Defaults to False.
        show_plot (bool, optional): Whether to show the plot. Defaults to True.
        output_dir (str, optional): Directory to save plot. Defaults to current directory.

    Returns:
        matplotlib.figure.Figure: The generated figure object
    """
    if Ns_all is None:
        Ns_all = np.arange(0, 9)

    plt.figure(figsize=(10, 6))

    for user in range(2, 5):  # users 2, 3, 4
        K = user - 1
        obj_mean = np.mean(np.squeeze(Obj_all[K, :, :]), axis=0)
        # Only plot if there's actually data (non-zero)
        if np.any(obj_mean != 0):
            plt.plot(Ns_all, obj_mean, 'o-', color=COLORS[K-1],
                     linewidth=2, markersize=6,
                     label=f'{user} users (K={K})')

    plt.xlabel('Number of Sensing Streams (Ns)')
    plt.ylabel('Objective Value')
    plt.title('Performance vs Number of Sensing Streams')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if save_plots:
        plot_file = os.path.join(output_dir, 'fig3_performance_vs_Ns_SCA.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {plot_file}")

    if show_plot:
        plt.show()
    else:
        plt.close()

    return plt.gcf()


def main():
    """Main function when run as a script."""
    parser = argparse.ArgumentParser(
        description='Plot performance vs number of sensing streams (Ns)'
    )
    parser.add_argument('--data-file', type=str, default='data.mat',
                        help='Path to data.mat file')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display the plot')
    parser.add_argument('--save', action='store_true',
                        help='Save the plot to file')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save output files')
    
    args = parser.parse_args()

    try:
        data = scipy.io.loadmat(args.data_file)
        Obj_all = data['Obj_all']
        plot_performance_vs_Ns(
            Obj_all,
            save_plots=args.save,
            show_plot=not args.no_show,
            output_dir=args.output_dir
        )
    except Exception as e:
        print(f"Error plotting data: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
