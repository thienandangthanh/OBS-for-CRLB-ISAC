#!/usr/bin/env python3
"""
Python translation of proposed_SCA.m for Figure 5 performance vs K analysis.

This script implements the Successive Convex Approximation (SCA) algorithm
for the OBS-for-CRLB-ISAC system, specifically for analyzing the performance
versus number of communication users K.

Original MATLAB file: Fig_5performance_vs_K/proposed_SCA.m
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

    # Iterative refinement (5 iterations as in MATLAB)
    for _ in range(5):
        FIM = calculateFIM(L, noise_s, Ws, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        CRBM = np.linalg.inv(FIM)
        Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        Ws = (Q + Q.conj().T) @ Ws
        Ws = Ws * np.sqrt(1 / np.trace(Ws @ Ws.conj().T))

    return Ws


def main():
    """Main function for Figure 5 performance vs K analysis."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='OBS-for-CRLB-ISAC Figure 5 Performance vs K Analysis',
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
    config = config.get_scenario_config('fig5')

    print("Figure 5 Performance vs K Analysis")
    print("=" * 50)
    print(config)
    print()

    # System parameters from MATLAB code
    M = config.M  # number of targets
    K_max = config.K_max  # maximum number of users

    # Antenna configuration
    Nth = config.Nth  # transmit horizontal antennas
    Ntv = config.Ntv  # transmit vertical antennas
    Nt = config.Nt    # total transmit antennas

    Nrh = config.Nrh  # receive horizontal antennas
    Nrv = config.Nrv  # receive vertical antennas
    Nr = config.Nr # total receive antennas

    # Power and noise configuration
    Pt = config.Pt  # transmit power in linear scale (dBm to W)
    noise_c = config.noise_c  # communication noise power
    noise_s = config.noise_s  # sensing noise power
    L = config.L  # number of sensing snapshots
    kappa = config.kappa

    # Algorithm parameters
    tolerance = config.tolerance
    max_iterations = config.max_iterations

    # Simulation parameters
    I_in = config.I_in  # inner loop iterations (K parameter sweep)
    I_out = config.I_out  # outer loop iterations (channel realizations)

    # Initialize result storage
    SR_all = np.zeros((I_out, I_in))
    CRB_all = np.zeros((I_out, I_in))
    Time_all = np.zeros((I_out, I_in))
    Obj_all = np.zeros((I_out, I_in))

    # Total number of parameter combinations
    total_combinations = I_out * I_in

    print(f'Starting parameter sweep over {total_combinations} combinations ({I_out} x {I_in})')
    print(f'K range: 2, 4, 6, 8, 10, 12 users')
    print()

    # Fixed trade-off parameters
    delta_s = config.delta_s
    delta_c = config.delta_c

    # Number of sensing streams equals number of transmit antennas
    num_sensing_streams = config.num_sensing_streams

    # Set random seed for reproducibility
    rng = np.random.default_rng(1)

    # Generate target parameters (consistent across all scenarios)
    alpha = config.generate_alpha(rng)

    theta, phi = config.generate_target_angles(rng)

    print(f"Target angles (theta): {theta}")
    print(f"Target angles (phi): {phi}")
    print(f"Target reflection coefficients (alpha): {alpha}")
    print()

    # Generate all channel realizations (I_out x Nt x K_max)
    H_all = 1/np.sqrt(2) * (rng.standard_normal((I_out, Nt, K_max)) +
        1j * rng.standard_normal((I_out, Nt, K_max)))

    # Construct steering matrices (same for all scenarios)
    A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
    B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)
    U = np.diag(alpha)

    # Store convergence data for the last scenario (for debugging)
    Con_last = []

    # Main parameter sweep loop
    for k_par in range(1, total_combinations + 1):

        # Display progress
        if (k_par - 1) % 50 == 0 or k_par == 1:
            print(f'Processing combination {k_par}/{total_combinations} ({100*k_par/total_combinations:.1f}%)')

        # Map linear index to 2D indices (following MATLAB indexing)
        channel = ((k_par - 1) % I_out) + 1
        weight = ((k_par - 1) // I_out) + 1

        print(f'[channel={channel}, weight={weight}]')

        # Set number of communication users K based on weight
        K = weight * 2  # K = 2, 4, 6, 8, 10, 12

        # Extract channel realization for K users
        H = np.squeeze(H_all[channel - 1, :, :K])  # Convert to 0-based indexing

        # Start timing for this scenario
        start_time = perf_counter()

        # Initialize beamforming matrices
        # Communication beamforming: delta_c*H./vecnorm(H)
        H_norms = np.linalg.norm(H, axis=0)
        Wc = delta_c * H / H_norms

        # Initialize sensing beamforming matrix with random values
        # Ws = np.random.randn(Nt, num_sensing_streams) + 1j * np.random.randn(Nt, num_sensing_streams)
        Ws = rng.standard_normal((Nt, num_sensing_streams)) + 1j * rng.standard_normal((Nt, num_sensing_streams))

        # Alternative: use the initial_Ws function (commented in MATLAB)
        # Ws = initial_Ws(L, noise_s, Nt, num_sensing_streams, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

        # Combine beamforming matrices
        W = np.hstack([Wc, Ws])

        # Normalize to satisfy power constraint
        W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

        # Calculate initial FIM
        FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        W_last = W.copy()
        Con = []

        # SCA algorithm main loop
        for count in range(max_iterations):

            # Calculate T_k (total received power for each user)
            H_W_comm = H.conj().T @ W[:, :K]  # Communication part only
            T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)

            # Calculate alpha_k and beta_k for linearization
            alpha_k = T_k / (T_k - square_abs(np.diag(H_W_comm))) - 1
            beta_k = np.sqrt(1 + alpha_k) * np.diag(H_W_comm) / T_k

            # Construct linearization matrices
            Sigma1 = np.diag(np.sqrt(1 + alpha_k) * beta_k)
            Sigma2 = np.diag(square_abs(beta_k))

            # Calculate CRLB matrix
            CRBM = np.linalg.inv(FIM)

            # Construct matrix Q for sensing optimization
            Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

            # Construct optimization matrices
            C1 = np.hstack([delta_c * H @ Sigma1, np.zeros((Nt, num_sensing_streams))])
            C2 = 0.5 * delta_s * (Q + Q.conj().T) - delta_c * H @ Sigma2 @ H.conj().T

            # Dominant eigenvalue calculation (Type B construction)
            mu = np.abs(eigs(H @ Sigma2 @ H.conj().T, k=1, which='LM', return_eigenvectors=False)[0])
            C2 = delta_c * mu * np.eye(Nt) + C2

            # SCA update step (single iteration as in MATLAB)
            W = C1 + C2 @ W

            # Normalize to satisfy power constraint
            W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

            # Update FIM with new W
            FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

            # Calculate current objective value and metrics
            T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
            H_W_comm = H.conj().T @ W[:, :K]

            rate_p = np.log(T_k / (T_k - square_abs(np.diag(H_W_comm))))
            sum_rate = np.sum(rate_p)
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

        # Store convergence data for last scenario
        if k_par == total_combinations:
            Con_last = Con

        # Calculate final metrics
        end_time = perf_counter()
        Rx = W @ W.conj().T
        FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        CRB = np.linalg.inv(FIM)

        # Final sum rate calculation
        T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
        H_W_comm = H.conj().T @ W[:, :K]
        SR = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H_W_comm)))))

        # Store results
        # Use 2D indexing [channel-1, weight-1] because Python requires explicit 2D coordinates for 2D arrays (unlike MATLAB's linear indexing)
        CRB_all[channel-1, weight-1] = np.real(np.trace(np.linalg.inv(FIM)))
        SR_all[channel-1, weight-1] = SR
        Time_all[channel-1, weight-1] = end_time - start_time
        Obj_all[channel-1, weight-1] = delta_c * SR - delta_s * np.real(np.trace(np.linalg.inv(FIM)))

    # Save results to file
    if args.save_data:
        os.makedirs(args.output_dir, exist_ok=True)
        output_file = os.path.join(args.output_dir, 'data_SCA_Ns=3M.mat')
        sio.savemat(output_file, {
            'SR_all': SR_all,
            'CRB_all': CRB_all,
            'Time_all': Time_all,
            'Obj_all': Obj_all,
            'Con_last': np.array(Con_last) if Con_last else np.array([]),
            'I_out': I_out,
            'I_in': I_in,
            'delta_s': delta_s,
            'delta_c': delta_c,
            'K_values': np.arange(2, 2*I_in+1, 2),  # [2, 4, 6, 8, 10, 12]
            'alpha': alpha,
            'theta': theta,
            'phi': phi,
            'Nt': Nt,
            'Nr': Nr,
            'M': M,
            'L': L,
            'Pt': Pt,
            'noise_c': noise_c,
            'noise_s': noise_s
        })
        print(f'\nResults saved to: {output_file}')

    # Display results summary
    print(f'\nSimulation completed in {np.sum(Time_all):.2f} seconds')
    print(f'Average processing time per scenario: {np.mean(Time_all):.4f} seconds')
    print()
    print('Results summary (averaged over channels):')
    print('K\tObjective\tSum Rate\tCRLB Trace')
    print('-' * 50)
    for i in range(I_in):
        K = (i + 1) * 2
        obj_mean = np.mean(Obj_all[:, i])
        sr_mean = np.mean(SR_all[:, i])
        crb_mean = np.mean(CRB_all[:, i])
        print(f'{K}\t{obj_mean:.4f}\t\t{sr_mean:.4f}\t\t{crb_mean:.4f}')

    # Optional: Create basic plots (similar to MATLAB figure commands)
    try:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Plot 1: Objective, Sum Rate, and CRB vs K
        K_values = np.arange(2, 2*I_in+1, 2)
        ax1.plot(K_values, np.mean(Obj_all, axis=0), 'o-', label='Objective')
        ax1.plot(K_values, np.mean(SR_all, axis=0), 's-', label='Sum Rate')
        ax1.plot(K_values, np.mean(CRB_all, axis=0), '^-', label='CRB Trace')
        ax1.set_xlabel('Number of Communication Users')
        ax1.set_ylabel('Value')
        ax1.set_title('Performance vs Number of Users')
        ax1.legend()
        ax1.grid(True)

        # Plot 2: Processing time vs K
        ax2.plot(K_values, np.mean(Time_all, axis=0), 'o-')
        ax2.set_xlabel('Number of Communication Users')
        ax2.set_ylabel('Processing Time (seconds)')
        ax2.set_title('Computational Complexity vs Number of Users')
        ax2.grid(True)

        plt.tight_layout()

        if args.save_plots:
            plt.savefig(os.path.join(args.output_dir, 'figure5_performance_vs_K.png'), dpi=300, bbox_inches='tight')
            print(f'Plots saved to: {os.path.join(args.output_dir, "figure5_performance_vs_K.png")}')
        else:
            plt.show()

    except ImportError:
        print('Matplotlib not available, skipping plots')

    print('\nAnalysis complete!')


if __name__ == "__main__":
    main()
