#!/usr/bin/env python3
"""
Python translation of proposed_SCA.m for Figure 2 trade-off region analysis.

This script implements the Successive Convex Approximation (SCA) algorithm
for the OBS-for-CRLB-ISAC system, specifically for analyzing the trade-off
region between communication sum rate and sensing performance (CRLB).

Original MATLAB file: Fig2_trade_off_region/proposed_SCA.m
Translated to Python with equivalent functionality.
"""

import numpy as np
import argparse
import scipy.io as sio
from time import perf_counter
from scipy.sparse.linalg import eigs
import sys
import os

# Add the parent directory to the path to import utils
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils import (
    construct_steer_matrix_and_derivative_steer_matrix,
    calculateFIM,
    construct_matrixQ,
    SimulationConfig,
    db2pow,
    square_abs
)


def initial_Ws(L, noise_s, Nt, num_sensing_streams, A, dAtheta, dAphi, B, dBtheta, dBphi, U):
    """
    Initialize sensing beamforming matrix using iterative algorithm.

    Args:
        L: Number of sensing snapshots
        noise_s: Sensing noise power
        Nt: Number of transmit antennas
        num_sensing_streams: Number of sensing streams
        A, dAtheta, dAphi: Transmit steering matrices and derivatives
        B, dBtheta, dBphi: Receive steering matrices and derivatives
        U: Target reflection coefficient matrix

    Returns:
        Ws: Initial sensing beamforming matrix
    """
    # Initialize with random complex matrix
    Ws = np.random.randn(Nt, num_sensing_streams) + 1j * np.random.randn(Nt, num_sensing_streams)

    # Iterative refinement
    for _ in range(5):
        FIM = calculateFIM(L, noise_s, Ws, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        CRBM = np.linalg.inv(FIM)
        Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        Ws = (Q + Q.conj().T) @ Ws
        Ws = Ws * np.sqrt(1 / np.trace(Ws @ Ws.conj().T))

    return Ws


def main():
    """Main function for Figure 2 trade-off region analysis."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='OBS-for-CRLB-ISAC Figure 2 Trade-off Region Analysis',
        parents=[SimulationConfig.create_argument_parser()],
        conflict_handler='resolve'
    )

    # Add figure-specific arguments
    parser.add_argument('--save-data', action='store_true', default=True,
                        help='Save results to data files')
    parser.add_argument('--save-plots', action='store_true',
                        help='Save plots as files instead of displaying')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save output files')

    args = parser.parse_args()

    # Initialize configuration from arguments
    config = SimulationConfig.from_args(args)
    config = config.get_scenario_config('fig2')

    # Print configuration
    print("Figure 2 Trade-off Region Analysis")
    print("=" * 50)
    print(config)
    print()

    # Extract configuration parameters
    # Antenna configuration
    Nth = config.Nth
    Ntv = config.Ntv
    Nt = config.Nt

    Nrh = config.Nrh
    Nrv = config.Nrv
    Nr = config.Nr

    # System configuration
    K = config.K
    M_max = config.M_max

    # Algorithm configuration
    num_sensing_streams = config.num_sensing_streams
    tolerance = config.tolerance
    max_iterations = config.max_iterations

    # Power and noise configuration
    Pt = config.Pt
    noise_c = config.noise_c
    noise_s = config.noise_s
    L = config.L
    kappa = config.kappa

    # Sweep configuration
    delta_all = config.get_delta_range()  # -7:0.2:4.8

    I_in = len(delta_all)
    I_out = config.I_out

    # Initialize result storage
    SR_all = np.zeros((I_out, I_in))
    CRB_all = np.zeros((I_out, I_in))
    Time_all = np.zeros((I_out, I_in))

    # Total number of parameter combinations
    total_combinations = I_out * I_in

    print(f'Starting parameter sweep over {total_combinations} combinations ({I_out} x {I_in})')
    print(f'Delta range: [{delta_all[0]:.1f}, {delta_all[-1]:.1f}] dB with {delta_all[1]-delta_all[0]:.1f} dB steps')

    # Fixed number of targets for this analysis
    m_target = 2
    M = m_target

    # Generate all random parameters once at the beginning for consistency
    rng = np.random.default_rng(1)
    alpha_all = config.alpha_base * (1 + config.alpha_variance * rng.standard_normal(M_max)) * \
        np.exp(1j * 2 * np.pi * rng.random(M_max))

    theta_range = config.theta_range
    phi_range = config.phi_range
    theta_all = theta_range[0] + (theta_range[1] - theta_range[0]) * rng.random(M_max)
    phi_all = phi_range[0] + (phi_range[1] - phi_range[0]) * rng.random(M_max)

    H_all = 1/np.sqrt(2) * (rng.standard_normal((I_out, Nt, K)) +
        1j * rng.standard_normal((I_out, Nt, K)))

    # Extract target parameters
    alpha = alpha_all[:m_target]
    theta = theta_all[:m_target]
    phi = phi_all[:m_target]

    print(f"Target angles (theta): {theta}")
    print(f"Target angles (phi): {phi}")
    print(f"Target reflection coefficients (alpha): {alpha}")
    print()

    # Construct steering matrices (same for all scenarios)
    A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
    B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)
    U = np.diag(alpha)

    # Store convergence data for the last scenario (for plotting)
    Con_last = []

    # Main parameter sweep loop
    for k_par in range(1, total_combinations + 1):
        
        # Display progress
        if (k_par - 1) % 100 == 0 or k_par == 1:
            print(f'Processing combination {k_par}/{total_combinations} ({100*k_par/total_combinations:.1f}%)')

        # Map linear index to 2D indices
        channel = ((k_par - 1) % I_out) + 1
        weight = ((k_par - 1) // I_out) + 1

        # Set trade-off parameters
        delta_s = config.delta_s
        delta_c = 10**(delta_all[weight - 1])  # Convert to 0-based indexing

        # Extract channel realization
        H = H_all[channel - 1, :, :]  # Convert to 0-based indexing

        # Start timing for this scenario
        start_time = perf_counter()

        # Initialize beamforming matrices with fresh random values for each scenario
        rng_scenario = np.random.default_rng(k_par + 1000)  # Different seed for each scenario
        Wc = rng_scenario.standard_normal((Nt, K)) + 1j * rng_scenario.standard_normal((Nt, K))
        Ws = rng_scenario.standard_normal((Nt, num_sensing_streams)) + 1j * rng_scenario.standard_normal((Nt, num_sensing_streams))

        # Combine and normalize beamforming matrices
        W = np.hstack([Wc, Ws])
        W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

        # Initial FIM calculation
        FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        W_last = W.copy()
        Con = []

        # SCA iteration loop
        for count in range(max_iterations):
            
            # Update auxiliary variables for communication
            T_k = np.sum(square_abs(H.conj().T @ W[:, :K]), axis=1) + noise_c * np.ones(K)
            alpha_k = T_k / (T_k - square_abs(np.diag(H.conj().T @ W[:, :K]))) - 1
            beta_k = np.sqrt(1 + alpha_k) * np.diag(H.conj().T @ W[:, :K]) / T_k

            Sigma1 = np.diag(np.sqrt(1 + alpha_k) * beta_k)
            Sigma2 = np.diag(square_abs(beta_k))

            # Update FIM and construct matrix Q for sensing
            CRBM = np.linalg.inv(FIM)
            Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

            # Construct matrices C1 and C2 for SCA update
            C1 = np.hstack([delta_c * H @ Sigma1, np.zeros((Nt, num_sensing_streams), dtype=complex)])
            C2 = 0.5 * delta_s * (Q + Q.conj().T) - delta_c * H @ Sigma2 @ H.conj().T

            # Dominant eigenvalue calculation (Type B construction)
            mu = np.abs(eigs(H @ Sigma2 @ H.conj().T, k=1, which='LM', return_eigenvectors=False)[0])
            C2 = delta_c * mu * np.eye(Nt) + C2

            # SCA update step
            W = C1 + C2 @ W
            W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

            # Update FIM
            FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

            # Calculate performance metrics for convergence tracking
            T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
            obj = delta_c * np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))) - \
                delta_s * np.trace(np.linalg.inv(FIM))

            # Store convergence data
            sum_rate = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W)))))
            trace_inv_fim = -np.trace(np.linalg.inv(FIM))
            Con.append([sum_rate, trace_inv_fim, obj])

            # Check convergence
            if np.linalg.norm(W - W_last) < tolerance:
                break
            else:
                W_last = W.copy()

        # Final calculations for this scenario
        T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
        SR = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W)))))
        CRB_trace = np.trace(np.linalg.inv(FIM))
        computation_time = perf_counter() - start_time

        # Store results in proper matrix locations
        SR_all[channel - 1, weight - 1] = SR  # Convert to 0-based indexing
        CRB_all[channel - 1, weight - 1] = CRB_trace
        Time_all[channel - 1, weight - 1] = computation_time

        # Store convergence data for the last scenario (for plotting)
        if k_par == total_combinations:
            Con_last = np.array(Con)

    print('Parameter sweep completed!')
    print(f'Average computation time per scenario: {np.mean(Time_all):.3f} seconds')
    print(f'Total computation time: {np.sum(Time_all)/60:.1f} minutes')

    # Save results to file
    if args.save_data:
        output_dir = args.output_dir

        # Prepare all variables for saving (matching MATLAB approach)
        save_dict = {
            'SR_all': SR_all,
            'CRB_all': CRB_all,
            'Time_all': Time_all,
            'delta_all': delta_all,
            'I_out': I_out,
            'I_in': I_in,
            'Nth': Nth,
            'Ntv': Ntv,
            'Nt': Nt,
            'Nrh': Nrh,
            'Nrv': Nrv,
            'Nr': Nr,
            'K': K,
            'M': M,
            'Pt': Pt,
            'noise_c': noise_c,
            'noise_s': noise_s,
            'L': L,
            'tolerance': tolerance
        }

        # Save to .mat file for MATLAB compatibility
        sio.savemat(os.path.join(output_dir, 'data_proposed_SCA.mat'), save_dict)
        print('Results saved to data_proposed_SCA.mat')

    # Display summary statistics
    print('\nSummary Statistics:')
    print(f'Sum Rate - Min: {np.min(SR_all):.4f}, Max: {np.max(SR_all):.4f}, Mean: {np.mean(SR_all):.4f}')
    print(f'CRB Trace - Min: {np.min(CRB_all):.4e}, Max: {np.max(CRB_all):.4e}, Mean: {np.mean(CRB_all):.4e}')

    # Plot convergence for the last scenario
    if len(Con_last) > 0:
        try:
            import matplotlib.pyplot as plt

            plt.figure(1)
            plt.plot(Con_last[:, 0], color='blue', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Sum Rate (nat/s/Hz)')
            plt.title('Sum Rate Convergence (Last Scenario)')
            plt.grid(True)

            plt.figure(2)
            plt.plot(Con_last[:, 1], color='red', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Trace of inverse FIM')
            plt.title('Sensing Performance Convergence (Last Scenario)')
            plt.grid(True)

            plt.figure(3)
            plt.plot(Con_last[:, 2], color='green', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Total Objective Value')
            plt.title('Objective Function Convergence (Last Scenario)')
            plt.grid(True)

            plt.tight_layout()

            if args.save_plots:
                # Save plots as files
                plot_names = [
                    'fig2_sum_rate_vs_iteration.png',
                    'fig2_inverse_fim_trace_vs_iteration.png',
                    'fig2_objective_value_vs_iteration.png'
                ]

                for fig_num, plot_name in enumerate(plot_names, 1):
                    plot_file = os.path.join(args.output_dir, plot_name)
                    plt.figure(fig_num)
                    plt.tight_layout()
                    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                    print(f"Plot {fig_num} saved to: {plot_file}")
                plt.close('all')
            else:
                plt.show()

        except ImportError:
            print("Matplotlib not available. Skipping plots.")

    return SR_all, CRB_all, Time_all, Con_last


if __name__ == "__main__":
    main()
