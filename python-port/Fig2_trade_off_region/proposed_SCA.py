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
    parser.add_argument('--save-data', action='store_true',
                        help='Save results to data files')
    parser.add_argument('--save-plots', action='store_true',
                        help='Save plots as files instead of displaying')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save output files')
    parser.add_argument('--m-target', type=int, default=2,
                        help='Number of targets')
    parser.add_argument('--k-par', type=int, default=2,
                        help='Parameter index for sweep')

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
    # NOTE:differences from figure 1
    # no M, use M_max instead
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

    # NOTE: new variables, not exist in Figure 1
    I_in = len(delta_all)
    I_out = config.I_out

    # Initialize result storage
    SR_all = np.zeros((I_out, I_in))
    CRB_all = np.zeros((I_out, I_in))
    Time_all = np.zeros((I_out, I_in))

    # Main computation loop
    m_target = args.m_target  # Default: 2
    k_par = args.k_par        # Default: 2

    print(f"Running analysis for m_target={m_target}, k_par={k_par}")
    print(f"Delta range: {delta_all[0]:.1f} to {delta_all[-1]:.1f} with {len(delta_all)} points")
    print()

    M = m_target
    channel = ((k_par - 1) % I_out) + 1  # Convert to 1-based indexing
    weight = ((k_par - 1) // I_out) + 1  # Convert to 1-based indexing

    print(f"Channel: {channel}, Weight: {weight}")

    delta_s = config.delta_s
    # NOTE:
    # delta_c is declared differently in figure 1
    # delta_c = lin[de]
    delta_c = 10**(delta_all[weight - 1])  # Convert to 0-based indexing

    print(f"delta_s: {delta_s}, delta_c: {delta_c:.6e}")
    print()

    # Initialize random number generator with fixed seed for reproducibility
    # TODO: use twister instead
    # numpy uses PCG64 as default
    # NOTE: figure 1 set seed with value of channel
    rng = np.random.default_rng(1)

    # Generate random parameters
    alpha_all = config.alpha_base * (1 + config.alpha_variance * rng.standard_normal(M_max)) * \
        np.exp(1j * 2 * np.pi * rng.random(M_max))

    # NOTE:
    # cannot use config.generate_target_angles
    # theta_all, phi_all = config.generate_target_angles(rng)
    # because this figure use rng.random(M_max)
    theta_range = config.theta_range
    phi_range = config.phi_range
    theta_all = theta_range[0] + (theta_range[1] - theta_range[0]) * rng.random(M_max)
    phi_all = phi_range[0] + (phi_range[1] - phi_range[0]) * rng.random(M_max)

    # Generate all channel realizations
    # NOTE: H_all is different from MATLAB because of random
    # NOTE: H_all is different from figure 1
    # in figure 1, I_out is not used
    H_all = 1/np.sqrt(2) * (rng.standard_normal((I_out, Nt, K)) +
        1j * rng.standard_normal((I_out, Nt, K)))
    # print(f"H_all: {H_all}")

    # NOTE: actually, only one value of H is used
    # Extract parameters for current scenario
    H = H_all[channel - 1, :, :]  # Convert to 0-based indexing
    # print(f"H: {H}")
    alpha = alpha_all[:m_target]
    theta = theta_all[:m_target]
    phi = phi_all[:m_target]

    print(f"Target angles (theta): {theta}")
    print(f"Target angles (phi): {phi}")
    print(f"Target reflection coefficients (alpha): {alpha}")
    print()

    # Start timing
    start_time = perf_counter()

    # Construct steering matrices
    A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
    B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)

    U = np.diag(alpha)

    # Initialize beamforming matrices
    # Communication beamforming: random initialization
    # TODO: Wc is different from MATLAB
    Wc = rng.standard_normal((Nt, K)) + 1j * rng.standard_normal((Nt, K))
    # print(f"Wc: {Wc}")

    # Sensing beamforming: random initialization
    # TODO: Ws is different from MATLAB
    Ws = rng.standard_normal((Nt, num_sensing_streams)) + 1j * rng.standard_normal((Nt, num_sensing_streams))
    # print(f"Ws: {Ws}")

    # Alternative initialization using initial_Ws (commented out as in MATLAB)
    # Ws = initial_Ws(L, noise_s, Nt, num_sensing_streams, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

    # Combine communication and sensing beamforming
    W = np.hstack([Wc, Ws])

    # Normalize to satisfy power constraint
    W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

    # Initial FIM calculation
    FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
    W_last = W.copy()
    # print(f'FIM: {FIM}')

    # Convergence tracking
    Con = []

    print("Starting SCA iterations...")

    # Main SCA iteration loop
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

        C2 = -0.5 * delta_s * (Q + Q.conj().T) + delta_c * H @ Sigma2 @ H.conj().T

        # Fig2 uses Type B C2 construction (different from Fig1)
        # Dominant eigenvalue calculation
        # MATLAB
        # mu=abs(eigs(H*Sigma2*H',1,'LM'));
        mu = np.abs(eigs(H @ Sigma2 @ H.conj().T, k=1, which='LM', return_eigenvectors=False)[0])
        C2 = delta_c * mu * np.eye(Nt) - C2

        # SCA update step
        W = C1 + C2 @ W

        # Normalize to satisfy power constraint
        W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

        # Update FIM
        FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

        # Calculate performance metrics
        T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
        rate_p = np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))

        # Total objective value
        obj = delta_c * np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))) - \
            delta_s * np.trace(np.linalg.inv(FIM))

        # Store convergence data
        sum_rate = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W)))))
        trace_inv_fim = -np.trace(np.linalg.inv(FIM))

        Con.append([sum_rate, trace_inv_fim, obj])

        # Check convergence
        if np.linalg.norm(W - W_last) < tolerance:
            print(f"Converged after {count + 1} iterations")
            break
        else:
            W_last = W.copy()

    else:
        print(f"Reached maximum iterations ({max_iterations})")

    # Convert convergence data to numpy array
    Con = np.array(Con)

    # Final calculations
    Rx = W @ W.conj().T
    FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
    CRB = np.linalg.inv(FIM)
    SR = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W)))))

    # Store final results - ensure real values for storage
    # CRB_trace = np.real(np.trace(np.linalg.inv(FIM)))  # Extract real part to avoid complex warning
    CRB_trace = np.trace(np.linalg.inv(FIM))
    CRB_all[k_par - 1] = CRB_trace
    SR_all[k_par - 1] = SR
    computation_time = perf_counter() - start_time
    Time_all[k_par - 1] = computation_time

    # Print final results
    print("\nFinal Results:")
    print(f"Sum Rate (SR): {SR:.6f} nat/s/Hz")
    print(f"Trace of inverse FIM (CRB): {CRB_trace:.6e}")  # Use scalar value for formatting
    print(f"Final objective value: {np.real(Con[-1, 2]):.6f}")  # Extract real part for formatting
    print(f"Computation time: {computation_time:.3f} seconds")  # Use scalar variable
    print()

    # Save results if requested
    if args.save_data:
        output_dir = args.output_dir

        # Save main results
        results_dict = {
            'SR_all': SR_all,
            'CRB_all': CRB_all,
            'Time_all': Time_all,
            'Con': Con,
            'config': {
                'delta_s': delta_s,
                'delta_c': delta_c,
                'delta_all': delta_all,
                'K': K,
                'M': M,
                'Nt': Nt,
                'Nr': Nr,
                'L': L,
                'tolerance': tolerance,
                'max_iterations': max_iterations
            }
        }

        # Save to .mat file for MATLAB compatibility
        sio.savemat(os.path.join(output_dir, 'data_proposed_SCA_python.mat'), results_dict)

        # Save to .npz file for Python
        # np.savez(os.path.join(output_dir, 'data_proposed_SCA_python.npz'), **results_dict)

        print(f"Results saved to {output_dir}/")
        print("Files created:")
        print("  - data_proposed_SCA_python.mat (MATLAB format)")
        # print("  - data_proposed_SCA_python.npz (Python format)")

    # Optional: Quick visualization of convergence
    try:
        import matplotlib.pyplot as plt

        plt.figure(1)
        plt.plot(Con[:, 0], color='blue')
        plt.xlabel('Iteration')
        plt.ylabel('Sum Rate (nat/s/Hz)')
        plt.title('Sum Rate Convergence')
        plt.grid(True)

        plt.figure(2)
        plt.plot(Con[:, 1], color='red')
        plt.xlabel('Iteration')
        plt.ylabel('Trace of inverse FIM')
        plt.title('Sensing Performance Convergence')
        plt.grid(True)

        plt.figure(3)
        plt.plot(Con[:, 2], color='green')
        plt.xlabel('Iteration')
        plt.ylabel('Total Objective Value')
        plt.title('Objective Function Convergence')
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

    return SR_all, CRB_all, Time_all, Con


if __name__ == "__main__":
    main()
