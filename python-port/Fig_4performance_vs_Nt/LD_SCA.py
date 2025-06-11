#!/usr/bin/env python3
"""
Python translation of LD_SCA.m for Figure 4 performance vs Nt analysis.

This script implements the LD (Low-Dimensional) SCA algorithm for the
OBS-for-CRLB-ISAC system, specifically for analyzing the performance
versus number of transmit antennas (Nt).

The LD-SCA uses an extended variable formulation with projection to
handle the optimization constraints more efficiently.

Original MATLAB file: Fig_4performance_vs_Nt/LD_SCA.m
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


def projection_ellipsoid(A, A_p, X0, Pt):
    """
    Project the extended variable X0 onto the power constraint ellipsoid.

    Args:
        A: Extended system matrix RS
        A_p: Pseudo-inverse of RS (RS_p)
        X0: Extended variable to project
        Pt: Total power constraint

    Returns:
        X: Projected extended variable
    """
    W = A @ X0
    W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))
    X = A_p @ W
    return X


def main():
    """Main function for Figure 4 LD-SCA performance vs Nt analysis."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='OBS-for-CRLB-ISAC Figure 4 LD-SCA Performance vs Nt Analysis',
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
    print("Figure 4 LD-SCA Performance vs Nt Analysis")
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

    print(f'Starting LD-SCA antenna scaling analysis: {I_out} realizations x {I_in} antenna configurations')
    print(f'Antenna configurations: Nt = [8, 16, 32, 64, 128]')
    print()

    # Main parameter sweep loop - for testing, just run first scenario (k_par=150 as in MATLAB)
    for k_par in [150]:  # Original runs k_par=150 only

        channel = ((k_par - 1) % I_out) + 1
        weight = ((k_par - 1) // I_out) + 1

        print(f'Processing: channel={channel}, weight={weight}')

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
        H = np.squeeze(H_all[channel-1, :, :])  # Extract specific channel realization

        # Start timing for this scenario
        start_time = perf_counter()

        # Construct steering matrices
        A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
        B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)

        U = np.diag(alpha)

        # Construct extended system matrix RS = [H, A, dAtheta, dAphi]
        RS = np.hstack([H, A, dAtheta, dAphi])
        RS_p = np.linalg.pinv(RS)
        RS_co = RS.conj().T @ RS

        # Extended channel and steering matrices
        H_e = RS.conj().T @ H
        A_e = RS.conj().T @ A
        dAtheta_e = RS.conj().T @ dAtheta
        dAphi_e = RS.conj().T @ dAphi

        # Initialize beamforming matrices
        # Communication beamforming: delta_c*H./vecnorm(H)
        H_norms = np.linalg.norm(H, axis=0)
        Wc = delta_c * H / H_norms

        # Random sensing beamforming initialization
        Ws = rng.standard_normal((Nt, num_sensing_streams)) + \
            1j * rng.standard_normal((Nt, num_sensing_streams))

        # Combine beamforming matrices
        W = np.hstack([Wc, delta_s * Ws])

        # Normalize total power
        W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

        # Convert to extended variable P
        P = RS_p @ W

        # Initial FIM calculation
        FIM = calculateFIM(L, noise_s, P @ P.conj().T, A_e, dAtheta_e, dAphi_e, B, dBtheta, dBphi, U)
        P_last = P.copy()
        Con = []

        # SCA iteration loop
        for count in range(max_iterations):

            # Update auxiliary variables for communication
            T_k = np.sum(square_abs(H_e.conj().T @ P[:, :K]), axis=1) + noise_c * np.ones(K)
            alpha_k = T_k / (T_k - square_abs(np.diag(H_e.conj().T @ P[:, :K]))) - 1
            beta_k = np.sqrt(1 + alpha_k) * np.diag(H_e.conj().T @ P[:, :K]) / T_k

            Sigma1 = np.diag(np.sqrt(1 + alpha_k) * beta_k)
            Sigma2 = np.diag(square_abs(beta_k))

            # Update FIM and construct matrix Q for sensing
            CRBM = np.linalg.inv(FIM)
            Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A_e, dAtheta_e, dAphi_e, B, dBtheta, dBphi, U)

            # Construct matrices C1 and C2 for SCA update
            C1 = np.hstack([delta_c * H_e @ Sigma1, np.zeros((RS.shape[1], num_sensing_streams), dtype=complex)])
            C2 = 0.5 * delta_s * (Q + Q.conj().T) - delta_c * H_e @ Sigma2 @ H_e.conj().T

            # Dominant eigenvalue calculation
            mu = np.abs(eigs(H_e @ Sigma2 @ H_e.conj().T, k=1, which='LM', return_eigenvectors=False)[0])
            C2 = delta_c * mu * RS_co + C2

            # Inner loop iterations (20 iterations as in MATLAB)
            for iter_inner in range(20):
                Linear = C1 + C2 @ P
                P = np.linalg.solve(RS_co, Linear)
                P = projection_ellipsoid(RS, RS_p, P, Pt)

            # Update FIM
            FIM = calculateFIM(L, noise_s, P @ P.conj().T, A_e, dAtheta_e, dAphi_e, B, dBtheta, dBphi, U)

            # Calculate performance metrics for convergence tracking
            # NOTE: why sum with axis=1?
            # MATLAB: T_k=sum(square_abs(H_e'*P),2)+noise_c*ones(K,1);
            T_k = np.sum(square_abs(H_e.conj().T @ P), axis=1) + noise_c * np.ones(K)
            sum_rate = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H_e.conj().T @ P)))))
            trace_inv_fim = np.real(np.trace(np.linalg.inv(FIM)))
            obj = delta_c * sum_rate - delta_s * trace_inv_fim

            # Store convergence data
            trace_inv_fim_neg = -trace_inv_fim
            Con.append([sum_rate, trace_inv_fim_neg, obj])

            # Check convergence
            if np.linalg.norm(RS @ P - RS @ P_last) < tolerance:
                break
            else:
                P_last = P.copy()

        # Final calculations for this scenario
        W = RS @ P
        FIM = calculateFIM(L, noise_s, W @ W.conj().T, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        CRB = np.linalg.inv(FIM)

        # Calculate sum rate for communication users only
        SR = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W[:, :K])))))

        computation_time = perf_counter() - start_time

        # Store results
        CRB_all[channel-1, weight-1] = np.real(np.trace(np.linalg.inv(FIM)))
        SR_all[channel-1, weight-1] = SR
        Time_all[channel-1, weight-1] = computation_time
        Obj_all[channel-1, weight-1] = delta_c * SR - delta_s * np.real(np.trace(np.linalg.inv(FIM)))

        print(f'Completed scenario {k_par}: SR={SR:.4f}, CRB_trace={CRB_all[channel-1, weight-1]:.4e}, Time={computation_time:.3f}s')

    print('\nLD-SCA analysis completed!')
    print(f'Total computation time: {np.sum(Time_all):.3f} seconds')

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
        file_name = 'data_SCA_LD_Ns=3M.mat'
        sio.savemat(os.path.join(output_dir, file_name), save_dict)
        print(f'Results saved to {file_name}')

    # Display summary statistics
    print('\nSummary Statistics:')
    valid_results = np.nonzero(SR_all)[0]
    if len(valid_results) > 0:
        print(f'Sum Rate - Min: {np.min(SR_all[valid_results]):.4f}, Max: {np.max(SR_all[valid_results]):.4f}, Mean: {np.mean(SR_all[valid_results]):.4f}')
        print(f'CRB Trace - Min: {np.min(CRB_all[valid_results]):.4e}, Max: {np.max(CRB_all[valid_results]):.4e}, Mean: {np.mean(CRB_all[valid_results]):.4e}')
        print(f'Objective - Min: {np.min(Obj_all[valid_results]):.4f}, Max: {np.max(Obj_all[valid_results]):.4f}, Mean: {np.mean(Obj_all[valid_results]):.4f}')
        print(f'Computation Time - Min: {np.min(Time_all[valid_results]):.3f}s, Max: {np.max(Time_all[valid_results]):.3f}s, Mean: {np.mean(Time_all[valid_results]):.3f}s')

    # Plot convergence for the scenario
    if len(Con) > 0:
        try:
            import matplotlib.pyplot as plt

            # Plot 1: Sum rate convergence
            plt.figure(1, figsize=(8, 6))
            Con = np.array(Con)
            plt.plot(Con[:, 0], color='blue', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Sum Rate (nat/s/Hz)')
            plt.title('LD-SCA Sum Rate Convergence')
            plt.grid(True)

            # Plot 2: Trace of inverse FIM convergence
            plt.figure(2, figsize=(8, 6))
            plt.plot(Con[:, 1], color='red', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Trace of inverse FIM')
            plt.title('LD-SCA Sensing Performance Convergence')
            plt.grid(True)

            # Plot 3: Objective function convergence
            plt.figure(3, figsize=(8, 6))
            plt.plot(Con[:, 2], color='green', linewidth=1.5)
            plt.xlabel('Iteration')
            plt.ylabel('Total Objective Value')
            plt.title('LD-SCA Objective Function Convergence')
            plt.grid(True)

            if args.save_plots:
                # Save plots as files
                plot_names = [
                    'fig4_LD_SCA_sum_rate_convergence.png',
                    'fig4_LD_SCA_fim_convergence.png',
                    'fig4_LD_SCA_objective_convergence.png'
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
