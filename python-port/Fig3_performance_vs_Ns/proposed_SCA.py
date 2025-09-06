#!/usr/bin/env python3
"""
Python translation of proposed_SCA.m for Figure 3 performance vs Ns analysis.

This script implements the Successive Convex Approximation (SCA) algorithm
for the OBS-for-CRLB-ISAC system, specifically for analyzing the performance
versus number of sensing streams (Ns).

Original MATLAB file: Fig3_performance_vs_Ns/proposed_SCA.m
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
    """Main function for Figure 3 performance vs Ns analysis."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='OBS-for-CRLB-ISAC Figure 3 Performance vs Ns Analysis',
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
    config = config.get_scenario_config('fig3')

    # Print configuration
    print("Figure 3 Performance vs Ns Analysis")
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

    # System configuration - Fig3 specific values
    K_max = 4  # From MATLAB: K_max=4
    M_max = 3  # From MATLAB: M_max=3
    m_target = 3  # From MATLAB: m_target=3

    # Algorithm configuration
    tolerance = config.tolerance
    max_iterations = config.max_iterations

    # Power and noise configuration
    Pt = config.Pt
    noise_c = config.noise_c
    noise_s = config.noise_s
    L = config.L
    kappa = config.kappa

    # Sensing streams configuration
    Ns_all = [0, 1, 2, 3, 4, 5, 6, 7, 8]  # From MATLAB: Ns_all=[0,1,2,3,4,5,6,7,8]

    I_in = len(Ns_all)
    I_out = config.I_out

    # Initialize result storage
    SR_all = np.zeros(I_out * I_in)
    CRB_all = np.zeros(I_out * I_in)
    Time_all = np.zeros(I_out * I_in)
    # 4 to accommodate users 2,3,4
    # leave all elements of index 0 with zero value
    # store from index 1
    # user = 2, 3, 4 => index = 1, 2, 3
    # => K = index
    Obj_all = np.zeros((4, I_out, I_in))

    print(f'Ns range: {Ns_all}')
    print(f'Users to test: [2, 3, 4]')
    print(f'Combinations per user: {I_out} x {I_in} = {I_out * I_in}')
    print()

    # Main user loop (2, 3, 4 users)
    for user in range(2, 5):

        print(f'\nProcessing user configuration: {user} users')
        K = user - 1
        print(f'K = {K} communication users')

        # Main parameter sweep loop
        for k_par in range(1, I_out * I_in + 1):

            channel = ((k_par - 1) % I_out) + 1
            weight = ((k_par - 1) // I_out) + 1

            print(f'[User {user}] Processing: channel={channel}, weight={weight}')

            # Set trade-off parameters (fixed for Figure 3)
            delta_s = 1
            delta_c = 0.25

            num_sensing_streams = weight - 1  # weight-1 from Ns_all indexing

            # Skip if num_sensing_streams is negative (shouldn't happen with Ns_all=[0,1,2,...])
            if num_sensing_streams < 0:
                continue

            # Generate random parameters (same seed for consistency)
            rng = np.random.default_rng(1)
            alpha_all = 0.1 * (1 + 0.2 * rng.standard_normal(M_max)) * \
                np.exp(1j * 2 * np.pi * rng.random(M_max))

            theta_all = -np.pi/3 + 2*np.pi/3 * rng.random(M_max)
            phi_all = -np.pi/3 + 2*np.pi/3 * rng.random(M_max)

            # Different H_all dimension order: (Nt, K_max, I_out)
            H_all = 1/np.sqrt(2) * (rng.standard_normal((Nt, K_max, I_out)) +
                1j * rng.standard_normal((Nt, K_max, I_out)))

            # Extract channel and target parameters
            if K > 0:
                H = H_all[:, :K, channel-1]  # Select K users and specific channel
            else:
                H = np.zeros((Nt, 0), dtype=complex)  # No communication users case
            alpha = alpha_all[:m_target]
            theta = theta_all[:m_target]
            phi = phi_all[:m_target]

            # Start timing for this scenario
            start_time = perf_counter()

            # Construct steering matrices
            A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
            B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)

            U = np.diag(alpha)

            # Initialize beamforming matrices
            # Different Wc initialization: delta_c*H./vecnorm(H)
            if K > 0:
                H_norms = np.linalg.norm(H, axis=0)
                Wc = delta_c * H / H_norms
            else:
                Wc = np.zeros((Nt, 0), dtype=complex)  # No communication streams

            # Random Ws initialization (only if num_sensing_streams > 0)
            if num_sensing_streams > 0:
                Ws = rng.standard_normal((Nt, num_sensing_streams)) + 1j * rng.standard_normal((Nt, num_sensing_streams))
            else:
                Ws = np.zeros((Nt, 0), dtype=complex)

            # Combine beamforming matrices
            if K > 0 and num_sensing_streams > 0:
                W = np.hstack([Wc, Ws])
            elif K > 0:
                W = Wc.copy()
            elif num_sensing_streams > 0:
                W = Ws.copy()
            else:
                # This shouldn't happen in practice (no communication and no sensing)
                W = rng.standard_normal((Nt, 1)) + 1j * rng.standard_normal((Nt, 1))

            # Normalize total power
            W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

            # Initial FIM calculation
            FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
            W_last = W.copy()
            Con = []

            # SCA iteration loop
            for count in range(max_iterations):

                # Update auxiliary variables for communication
                if K > 0:
                    T_k = np.sum(square_abs(H.conj().T @ W[:, :K]), axis=1) + noise_c * np.ones(K)
                    alpha_k = T_k / (T_k - square_abs(np.diag(H.conj().T @ W[:, :K]))) - 1
                    beta_k = np.sqrt(1 + alpha_k) * np.diag(H.conj().T @ W[:, :K]) / T_k

                    Sigma1 = np.diag(np.sqrt(1 + alpha_k) * beta_k)
                    Sigma2 = np.diag(square_abs(beta_k))
                else:
                    # No communication users case
                    T_k = np.array([])
                    Sigma1 = np.zeros((0, 0), dtype=complex)
                    Sigma2 = np.zeros((0, 0), dtype=complex)

                # Update FIM and construct matrix Q for sensing
                CRBM = np.linalg.inv(FIM)
                Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

                # Construct matrices C1 and C2 for SCA update
                if K > 0 and num_sensing_streams > 0:
                    C1 = np.hstack([delta_c * H @ Sigma1, np.zeros((Nt, num_sensing_streams), dtype=complex)])
                elif K > 0:
                    C1 = delta_c * H @ Sigma1
                elif num_sensing_streams > 0:
                    C1 = np.zeros((Nt, num_sensing_streams), dtype=complex)
                else:
                    C1 = np.zeros((Nt, 1), dtype=complex)

                if K > 0:
                    C2 = 0.5 * delta_s * (Q + Q.conj().T) - delta_c * H @ Sigma2 @ H.conj().T
                else:
                    C2 = 0.5 * delta_s * (Q + Q.conj().T)

                # Dominant eigenvalue calculation (Type B construction)
                if K > 0:
                    mu = np.abs(eigs(H @ Sigma2 @ H.conj().T, k=1, which='LM', return_eigenvectors=False)[0])
                    C2 = delta_c * mu * np.eye(Nt) + C2
                else:
                    # No need for dominant eigenvalue in sensing-only case
                    pass

                # SCA update step (single iteration as in MATLAB)
                W = C1 + C2 @ W
                W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

                # Update FIM
                FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

                # Calculate performance metrics for convergence tracking
                if K > 0:
                    T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
                    sum_rate = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W)))))
                else:
                    sum_rate = 0.0

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

            # Calculate sum rate with proper indexing for communication users only
            if K > 0:
                SR = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W[:, :K])))))
            else:
                SR = 0.0

            computation_time = perf_counter() - start_time

            # Store results
            CRB_all[k_par-1] = np.real(np.trace(np.linalg.inv(FIM)))
            SR_all[k_par-1] = SR
            Time_all[k_par-1] = computation_time
            Obj_all[K, channel-1, weight-1] = delta_c * SR - delta_s * np.real(np.trace(np.linalg.inv(FIM)))

        # Save intermediate data after each user
        save_dict = {'Obj_all': Obj_all}
        sio.savemat(os.path.join(args.output_dir, 'data.mat'), save_dict)

    print('\nParameter sweep completed!')
    print(f'Average computation time per scenario: {np.mean(Time_all):.3f} seconds')
    print(f'Total computation time: {np.sum(Time_all)/60:.1f} minutes')

    # Save final results to file
    if args.save_data:
        output_dir = args.output_dir

        # Prepare all variables for saving (matching MATLAB approach)
        save_dict = {
            'Obj_all': Obj_all,
            'SR_all': SR_all,
            'CRB_all': CRB_all,
            'Time_all': Time_all,
            'Ns_all': Ns_all,
            'I_out': I_out,
            'I_in': I_in,
            'Nth': Nth,
            'Ntv': Ntv,
            'Nt': Nt,
            'Nrh': Nrh,
            'Nrv': Nrv,
            'Nr': Nr,
            'K_max': K_max,
            'M_max': M_max,
            'm_target': m_target,
            'Pt': Pt,
            'noise_c': noise_c,
            'noise_s': noise_s,
            'L': L,
            'tolerance': tolerance
        }

        # Save to .mat file for MATLAB compatibility
        sio.savemat(os.path.join(output_dir, 'data_proposed_SCA_fig3.mat'), save_dict)
        print('Results saved to data_proposed_SCA_fig3.mat')

    # Display summary statistics
    print('\nSummary Statistics:')
    print(f'Objective values - Min: {np.min(Obj_all):.4f}, Max: {np.max(Obj_all):.4f}, Mean: {np.mean(Obj_all):.4f}')
    print(f'Sum Rate - Min: {np.min(SR_all):.4f}, Max: {np.max(SR_all):.4f}, Mean: {np.mean(SR_all):.4f}')
    print(f'CRB Trace - Min: {np.min(CRB_all):.4e}, Max: {np.max(CRB_all):.4e}, Mean: {np.mean(CRB_all):.4e}')

    # Plot objective values vs Ns for different users
    try:
        from plot_SCA_performance_vs_Ns import plot_performance_vs_Ns
        plot_performance_vs_Ns(
            Obj_all,
            Ns_all=Ns_all,
            save_plots=args.save_plots,
            show_plot=not args.save_plots,
            output_dir=args.output_dir
        )
    except ImportError:
        print("Could not import plotting module. Skipping plots.")

    return Obj_all, SR_all, CRB_all, Time_all


if __name__ == "__main__":
    main() 
