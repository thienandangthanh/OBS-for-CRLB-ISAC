#!/usr/bin/env python3
"""
Python translation of proposed_SCA.m for Figure 4 performance vs Nt analysis.

This script implements the Successive Convex Approximation (SCA) algorithm
for the OBS-for-CRLB-ISAC system, specifically for analyzing the performance
versus number of transmit antennas (Nt).

Original MATLAB file: Fig_4performance_vs_Nt/proposed_SCA.m
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
    """Main function for Figure 4 performance vs Nt analysis."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='OBS-for-CRLB-ISAC Figure 4 Performance vs Nt Analysis',
        parents=[SimulationConfig.create_argument_parser()],
        conflict_handler='resolve'
    )

    # Add figure-specific arguments
    parser.add_argument('--save-data', action='store_true',
                        help='Save results to data files')
    parser.add_argument('--save-plots', action='store_true',
                        help='Save plots as files instead of displaying')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save output files')

    args = parser.parse_args()

    # Initialize configuration from arguments
    config = SimulationConfig.from_args(args)
    config = config.get_scenario_config('fig4')

    # Print configuration
    print("Figure 4 Performance vs Nt Analysis")
    print("=" * 50)
    print(config)
    print()

    # Extract configuration parameters
    # System configuration - Fig4 specific values
    K = config.K
    M = config.M

    # Algorithm configuration
    tolerance = config.tolerance
    max_iterations = config.max_iterations

    # Power and noise configuration
    Pt = config.Pt
    noise_c = config.noise_c
    noise_s = config.noise_s
    L = config.L
    kappa = config.kappa

    # Set trade-off parameters (fixed for Figure 4)
    delta_s = config.delta_s
    delta_c = config.delta_c

    # Sweep configuration
    I_in = config.I_in  # Number of antenna configurations
    I_out = config.I_out  # Number of channel realizations

    # Initialize result storage
    SR_all = np.zeros((I_out, I_in))
    CRB_all = np.zeros((I_out, I_in))
    Time_all = np.zeros((I_out, I_in))
    Obj_all = np.zeros((I_out, I_in))

    print(f'Starting antenna scaling analysis: {I_out} realizations x {I_in} antenna configurations')
    print(f'Antenna configurations: Nt = [8, 16, 32, 64, 128]')
    print()

    # Main parameter sweep loop
    for k_par in range(1, I_out * I_in + 1):

        channel = ((k_par - 1) % I_out) + 1
        weight = ((k_par - 1) // I_out) + 1

        print(f'Processing: channel={channel}, weight={weight} (Antenna config {weight}/5)')

        # Variable antenna configuration based on weight
        # weight = 1,2,3,4,5 corresponds to Nt = 8,16,32,64,128
        Nth = 4
        Ntv = 2**(weight + 1)
        Nt = Nth * Ntv

        # Receiver antenna configuration (scales with transmit)
        Nrh = 5
        Nrv = 2**(weight + 1)
        Nr = Nrh * Nrv

        # Fixed sensing streams
        num_sensing_streams = 3 * M  # 6 sensing streams

        # Generate random parameters (same seed for consistency)
        rng = np.random.default_rng(1)
        alpha = config.generate_alpha(rng)

        theta, phi = config.generate_target_angles(rng)

        # Channel matrix generation: (I_out, Nt, K)
        H_all = 1/np.sqrt(2) * (rng.standard_normal((I_out, Nt, K)) +
            1j * rng.standard_normal((I_out, Nt, K)))

        # Extract channel realization
        H = H_all[channel-1, :, :]  # Extract specific channel realization

        # Start timing for this scenario
        start_time = perf_counter()

        # Construct steering matrices
        A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
        B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)

        U = np.diag(alpha)

        # Initialize beamforming matrices
        # Communication beamforming: delta_c*H./vecnorm(H)
        H_norms = np.linalg.norm(H, axis=0)
        Wc = delta_c * H / H_norms

        # Random sensing beamforming initialization
        rng_scenario = np.random.default_rng(k_par + 1000)
        Ws = rng_scenario.standard_normal((Nt, num_sensing_streams)) + \
            1j * rng_scenario.standard_normal((Nt, num_sensing_streams))

        # Combine beamforming matrices
        W = np.hstack([Wc, Ws])

        # Normalize total power
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

            # SCA update step (single iteration as in MATLAB)
            W = C1 + C2 @ W
            W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

            # Update FIM
            FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

            # Calculate performance metrics for convergence tracking
            T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
            sum_rate = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W)))))
            trace_inv_fim = np.real(np.trace(np.linalg.inv(FIM)))
            obj = delta_c * sum_rate - delta_s * trace_inv_fim

            # Store convergence data
            trace_inv_fim_neg = -trace_inv_fim
            Con.append([sum_rate, trace_inv_fim_neg, obj])

            # Check convergence
            if np.linalg.norm(W - W_last) < tolerance:
                break
            else:
                W_last = W.copy()

        # Final calculations for this scenario
        Rx = W @ W.conj().T
        FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        CRB = np.linalg.inv(FIM)

        # Calculate sum rate for communication users only
        SR = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W[:, :K])))))

        computation_time = perf_counter() - start_time

        # Store results
        CRB_all[channel-1, weight-1] = np.real(np.trace(np.linalg.inv(FIM)))
        SR_all[channel-1, weight-1] = SR
        Time_all[channel-1, weight-1] = computation_time
        Obj_all[channel-1, weight-1] = delta_c * SR - delta_s * np.real(np.trace(np.linalg.inv(FIM)))

        # Display progress
        if (k_par-1) % 20 == 0 or k_par == I_out * I_in:
            progress = 100 * k_par / (I_out * I_in)
            avg_time = np.mean(Time_all[Time_all > 0])
            print(f'Progress: {progress:.1f}% - Avg time per scenario: {avg_time:.3f}s')

    print('\nParameter sweep completed!')
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
            'Obj_all': Obj_all,
            'I_out': I_out,
            'I_in': I_in,
            'K': K,
            'M': M,
            'Pt': Pt,
            'noise_c': noise_c,
            'noise_s': noise_s,
            'L': L,
            'tolerance': tolerance,
            'delta_s': delta_s,
            'delta_c': delta_c,
            'num_sensing_streams': num_sensing_streams
        }

        # Save to .mat file for MATLAB compatibility
        sio.savemat(os.path.join(output_dir, 'data_SCA_Ns=3M.mat'), save_dict)
        print('Results saved to data_SCA_Ns=3M.mat')

    # Display summary statistics
    print('\nSummary Statistics:')
    print(f'Sum Rate - Min: {np.min(SR_all):.4f}, Max: {np.max(SR_all):.4f}, Mean: {np.mean(SR_all):.4f}')
    print(f'CRB Trace - Min: {np.min(CRB_all):.4e}, Max: {np.max(CRB_all):.4e}, Mean: {np.mean(CRB_all):.4e}')
    print(f'Objective - Min: {np.min(Obj_all):.4f}, Max: {np.max(Obj_all):.4f}, Mean: {np.mean(Obj_all):.4f}')
    print(f'Computation Time - Min: {np.min(Time_all):.3f}s, Max: {np.max(Time_all):.3f}s, Mean: {np.mean(Time_all):.3f}s')

    # Plot convergence for the last scenario
    if len(Con) > 0:
        try:
            import matplotlib.pyplot as plt

            # Antenna sizes for plotting
            antenna_sizes = [8, 16, 32, 64, 128]

            # Plot 1: Sum rate convergence (last scenario)
            plt.figure(1, figsize=(8, 6))
            Con = np.array(Con)
            plt.plot(Con[:, 0], color='blue', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Sum Rate (nat/s/Hz)')
            plt.title('Sum Rate Convergence (Last Scenario)')
            plt.grid(True)

            # Plot 2: Trace of inverse FIM convergence (last scenario)
            plt.figure(2, figsize=(8, 6))
            plt.plot(Con[:, 1], color='red', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Trace of inverse FIM')
            plt.title('Sensing Performance Convergence (Last Scenario)')
            plt.grid(True)

            # Plot 3: Objective function convergence (last scenario)
            plt.figure(3, figsize=(8, 6))
            plt.plot(Con[:, 2], color='green', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Total Objective Value')
            plt.title('Objective Function Convergence (Last Scenario)')
            plt.grid(True)

            # Plot 4: Performance vs number of antennas
            plt.figure(4, figsize=(10, 6))
            
            # Calculate mean performance for each antenna configuration
            mean_obj = np.mean(Obj_all, axis=0)
            mean_sr = np.mean(SR_all, axis=0)
            mean_crb = np.mean(CRB_all, axis=0)
            mean_time = np.mean(Time_all, axis=0)

            plt.subplot(2, 2, 1)
            plt.plot(antenna_sizes, mean_obj, 'o-', linewidth=2)
            plt.xlabel('Number of Transmit Antennas')
            plt.ylabel('Objective Value')
            plt.title('Objective Function vs Nt')
            plt.grid(True)

            plt.subplot(2, 2, 2)
            plt.plot(antenna_sizes, mean_sr, 's-', linewidth=2, color='green')
            plt.xlabel('Number of Transmit Antennas')
            plt.ylabel('Sum Rate (nat/s/Hz)')
            plt.title('Sum Rate vs Nt')
            plt.grid(True)

            plt.subplot(2, 2, 3)
            plt.semilogy(antenna_sizes, mean_crb, '^-', linewidth=2, color='red')
            plt.xlabel('Number of Transmit Antennas')
            plt.ylabel('CRB Trace')
            plt.title('CRB Trace vs Nt')
            plt.grid(True)

            plt.subplot(2, 2, 4)
            plt.loglog(antenna_sizes, mean_time, 'd-', linewidth=2, color='purple')
            plt.xlabel('Number of Transmit Antennas')
            plt.ylabel('Average CPU Time (sec)')
            plt.title('Computation Time vs Nt')
            plt.grid(True)

            plt.tight_layout()

            if args.save_plots:
                # Save plots as files
                plot_names = [
                    'fig4_sum_rate_convergence.png',
                    'fig4_fim_convergence.png',
                    'fig4_objective_convergence.png',
                    'fig4_performance_vs_Nt.png'
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

    return SR_all, CRB_all, Time_all, Obj_all, Con


if __name__ == "__main__":
    main() 
