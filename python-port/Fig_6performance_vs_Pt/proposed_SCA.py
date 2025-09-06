#!/usr/bin/env python3
"""
Python translation of proposed_SCA.m for Figure 6 performance vs Pt analysis.

This script implements the Successive Convex Approximation (SCA) algorithm
for the OBS-for-CRLB-ISAC system, specifically for analyzing the performance
versus transmit power (Pt).

Original MATLAB file: Fig_6performance_vs_Pt/proposed_SCA.m
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
    """Main function for Figure 6 performance vs Pt analysis."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='OBS-for-CRLB-ISAC Figure 6 Performance vs Pt Analysis',
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
    config = config.get_scenario_config('fig6')

    # Print configuration
    print("Figure 6 Performance vs Pt Analysis")
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

    # System configuration - Fig6 specific values
    M = config.M  # Number of targets
    K = config.K  # Number of communication users

    # Algorithm configuration
    tolerance = config.tolerance
    max_iterations = config.max_iterations

    # Power and noise configuration (Pt varies, others are fixed)
    noise_c = config.noise_c
    noise_s = config.noise_s
    L = config.L
    # NOTE: kappa is not used
    kappa = config.kappa

    # Transmit power configuration
    I_in = config.I_in   # Number of power levels
    I_out = config.I_out  # Number of channel realizations

    # Initialize result storage
    SR_all = np.zeros((I_out, I_in))
    CRB_all = np.zeros((I_out, I_in))
    Time_all = np.zeros((I_out, I_in))
    Obj_all = np.zeros((I_out, I_in))

    print(f'Power levels: {I_in}')
    print(f'Channel realizations: {I_out}')
    print(f'Total combinations: {I_out * I_in}')
    print()

    # Set trade-off parameters (fixed for Figure 6)
    delta_s = config.delta_s   # No sensing weight - pure communication focus
    delta_c = config.delta_c   # Communication weight

    # Full antenna array for sensing
    num_sensing_streams = Nt

    # Main parameter sweep loop
    for k_par in range(1, I_out * I_in + 1):

        channel = ((k_par - 1) % I_out) + 1
        weight = ((k_par - 1) // I_out) + 1

        print(f'Processing: channel={channel}, weight={weight}')

        # Variable transmit power based on weight (Fig6 specific)
        Pt = db2pow(5 * weight - 10 - 30)  # dBm
        print(f'  Pt = {5 * weight - 10} dBm = {Pt:.2e} W')

        # Generate random parameters with fixed seed for reproducibility
        rng = np.random.default_rng(1)

        # Target parameters
        alpha = config.generate_alpha(rng)
        theta, phi = config.generate_target_angles(rng)

        # Channel matrix (different dimension order than Fig3: (I_out, Nt, K))
        H_all = 1/np.sqrt(2) * (rng.standard_normal((I_out, Nt, K)) +
            1j * rng.standard_normal((I_out, Nt, K)))

        # Extract specific channel realization (matching MATLAB squeeze(H_all(channel,:,1:K)))
        H = np.squeeze(H_all[channel-1, :, :K])  # Shape should be (Nt, K)

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

        # Sensing beamforming: random initialization
        Ws = rng.standard_normal((Nt, num_sensing_streams)) + 1j * rng.standard_normal((Nt, num_sensing_streams))

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
                print(f'  Converged at iteration {count+1}')
                break
            else:
                W_last = W.copy()

        # Final calculations for this scenario
        Rx = W @ W.conj().T
        FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

        # Calculate final sum rate (using only communication beamforming)
        T_k_final = np.sum(square_abs(H.conj().T @ W[:, :K]), axis=1) + noise_c * np.ones(K)
        SR = np.sum(np.log(T_k_final / (T_k_final - square_abs(np.diag(H.conj().T @ W[:, :K])))))

        computation_time = perf_counter() - start_time

        # Store results (using MATLAB indexing: k_par-1 maps to (channel-1, weight-1))
        CRB_all[channel-1, weight-1] = np.real(np.trace(np.linalg.inv(FIM)))
        SR_all[channel-1, weight-1] = SR
        Time_all[channel-1, weight-1] = computation_time
        Obj_all[channel-1, weight-1] = delta_c * SR - delta_s * np.real(np.trace(np.linalg.inv(FIM)))

        print(f'  Results: SR={SR:.4f}, CRB_trace={CRB_all[channel-1, weight-1]:.4e}, '
              f'Obj={Obj_all[channel-1, weight-1]:.4f}, Time={computation_time:.3f}s')

    print('\nParameter sweep completed!')
    print(f'Average computation time per scenario: {np.mean(Time_all):.3f} seconds')
    print(f'Total computation time: {np.sum(Time_all)/60:.1f} minutes')

    # Calculate and display mean performance across channels
    mean_obj = np.mean(Obj_all, axis=0)
    mean_sr = np.mean(SR_all, axis=0)
    mean_crb = np.mean(CRB_all, axis=0)

    print('\nMean Performance vs Transmit Power:')
    print('Weight | Pt (dBm) | Objective | Sum Rate | CRB Trace')
    print('-' * 55)
    for i in range(I_in):
        pt_dbm = 5 * (i + 1) - 10
        print(f'{i+1:6d} | {pt_dbm:8.1f} | {mean_obj[i]:9.4f} | {mean_sr[i]:8.4f} | {mean_crb[i]:9.4e}')

    # Save final results to file
    if args.save_data:
        output_dir = args.output_dir

        # Prepare all variables for saving (matching MATLAB approach)
        save_dict = {
            'SR_all': SR_all,
            'CRB_all': CRB_all,
            'Time_all': Time_all,
            'Obj_all': Obj_all,
            'I_in': I_in,
            'I_out': I_out,
            'M': M,
            'K': K,
            'Nth': Nth,
            'Ntv': Ntv,
            'Nt': Nt,
            'Nrh': Nrh,
            'Nrv': Nrv,
            'Nr': Nr,
            'noise_c': noise_c,
            'noise_s': noise_s,
            'L': L,
            'tolerance': tolerance,
            'delta_s': delta_s,
            'delta_c': delta_c,
            'num_sensing_streams': num_sensing_streams
        }

        # Save to .mat file for MATLAB compatibility
        file_name = 'data_SCA_Ns=3M.mat'
        sio.savemat(os.path.join(output_dir, file_name), save_dict)
        print(f'\nResults saved to {os.path.join(output_dir, file_name)}')

    # Display summary statistics
    print('\nSummary Statistics:')
    print(f'Objective values - Min: {np.min(Obj_all):.4f}, Max: {np.max(Obj_all):.4f}, Mean: {np.mean(Obj_all):.4f}')
    print(f'Sum Rate - Min: {np.min(SR_all):.4f}, Max: {np.max(SR_all):.4f}, Mean: {np.mean(SR_all):.4f}')
    print(f'CRB Trace - Min: {np.min(CRB_all):.4e}, Max: {np.max(CRB_all):.4e}, Mean: {np.mean(CRB_all):.4e}')

    # Plot performance vs Pt
    try:
        from plot_SCA_performance_vs_Pt import plot_performance_vs_Pt
        plot_performance_vs_Pt(
            Obj_all, SR_all, CRB_all,
            I_in=I_in,
            save_plots=args.save_plots,
            show_plot=not args.save_plots,
            output_dir=args.output_dir
        )
    except ImportError:
        print("Could not import plotting module. Skipping plots.")

    return Obj_all, SR_all, CRB_all, Time_all


if __name__ == "__main__":
    main()
