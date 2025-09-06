import numpy as np
from scipy.io import savemat
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt
from time import perf_counter
import os
import sys
import argparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.steering_matrix import construct_steer_matrix_and_derivative_steer_matrix
from utils.calculate_fim import calculateFIM
from utils.construct_matrixQ import construct_matrixQ
from utils.square_abs import square_abs
from utils.simulation_config import SimulationConfig

# Not used
def initial_Ws(L, noise_s, Nt, num_sensing_streams, A, dAtheta, dAphi, B, dBtheta, dBphi, U):
    """Initialize sensing beamforming matrix"""
    Ws = np.random.randn(Nt, num_sensing_streams) + 1j * np.random.randn(Nt, num_sensing_streams)

    for _ in range(5):
        FIM = calculateFIM(L, noise_s, Ws, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        CRBM = np.linalg.inv(FIM)
        Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        Ws = (Q + Q.conj().T) @ Ws
        Ws = Ws * np.sqrt(1 / np.trace(Ws @ Ws.conj().T))

    return Ws

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='OBS-for-CRLB-ISAC Figure 1 Convergence Analysis',
        parents=[SimulationConfig.create_argument_parser()],
        conflict_handler='resolve'
    )

    # Add figure-specific arguments
    parser.add_argument('--save-plots', action='store_true',
                        help='Save plots as files instead of displaying')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory to save output files')

    args = parser.parse_args()

    # Initialize configuration from arguments
    config = SimulationConfig.from_args(args)
    config = config.get_scenario_config('fig1')

    # Print configuration
    print(config)
    print()

    # Get configuration parameters
    Nth = config.Nth
    Ntv = config.Ntv
    Nt = config.Nt

    Nrh = config.Nrh
    Nrv = config.Nrv
    # NOTE:
    # Nr is not used
    # Nr = config.Nr

    # NOTE:
    # channel_number is not used
    # channel_number = config.channel_number
    K = config.K
    # NOTE:
    # No need to set M here
    # Because it's used internally in config.generate_target_angles
    # M = config.M

    num_sensing_streams = config.num_sensing_streams
    tolerance = config.tolerance
    max_iterations = config.max_iterations

    delta_s = config.delta_s
    lin = config.convergence_lin_values  # [0.05, 0.1, 0.15]

    Pt = config.Pt
    noise_c = config.noise_c
    noise_s = config.noise_s

    L = config.L
    # NOTE:
    # kappa is not used
    # kappa = config.kappa

    # Initialize random number generator
    rng = np.random.default_rng(config.random_seed)

    # Generate target angles using config method
    theta, phi = config.generate_target_angles(rng)

    delta_all = config.get_delta_range()
    # NOTE:
    # Per_all is not used
    Per_all = np.zeros((2, len(delta_all)))
    Con = [[], [], []]

    for channel in [2]:
        for de in range(3):
            delta_c = lin[de]

            # Generate alpha using config method
            rng = np.random.default_rng(channel)
            alpha = config.generate_alpha(rng)

            H = 1/np.sqrt(2) * (rng.standard_normal((Nt, K)) + 1j * rng.standard_normal((Nt, K)))

            # MATLAB uses tic
            start_time = perf_counter()

            A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
            B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)

            U = np.diag(alpha)

            Wc = rng.standard_normal((Nt, K)) + 1j * rng.standard_normal((Nt, K))
            Ws = rng.standard_normal((Nt, num_sensing_streams)) + 1j * rng.standard_normal((Nt, num_sensing_streams))

            # TODO: Why hstack is used?
            W = np.hstack([Wc, Ws])
            W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

            FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
            W_last = W.copy()

            for count in range(max_iterations):
                T_k = np.sum(square_abs(H.conj().T @ W[:, :K]), axis=1) + noise_c * np.ones(K)
                alpha_k = T_k / (T_k - square_abs(np.diag(H.conj().T @ W[:, :K]))) - 1
                beta_k = np.sqrt(1 + alpha_k) * np.diag(H.conj().T @ W[:, :K]) / T_k

                Sigma1 = np.diag(np.sqrt(1 + alpha_k) * beta_k)
                Sigma2 = np.diag(square_abs(beta_k))

                CRBM = np.linalg.inv(FIM)
                Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

                C1 = np.hstack([delta_c * H @ Sigma1, np.zeros((Nt, num_sensing_streams), dtype=complex)])
                C2 = -0.5 * delta_s * (Q + Q.conj().T) + delta_c * H @ Sigma2 @ H.conj().T

                # Use largest magnitude eigenvalue to match MATLAB's eigs(C2,1,'LM')
                # eigenvals = linalg.eigvals(C2)
                # mu = np.abs(eigenvals[np.argmax(np.abs(eigenvals))])
                mu = np.abs(eigs(C2, k=1, which='LM', return_eigenvectors=False)[0])
                C2 = mu * np.eye(Nt) - C2

                W = C1 + C2 @ W
                W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

                FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
                print(f'FIM: {FIM}')

                # TODO: Why axis=1
                T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
                # TODO: Why the dot of the first T_k. is omitted
                rate_p = np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))
                # TODO: Why the dot of the first T_k. is omitted
                obj = delta_c * np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))) - \
                    delta_s * np.trace(np.linalg.inv(FIM))

                # Append data to Con[de] to match MATLAB behavior
                Con[de].append([
                    # TODO: Why the dot of the first T_k. is omitted
                    np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))),
                    -np.trace(np.linalg.inv(FIM)),
                    obj
                ])

                # Check convergence
                norm_diff = np.linalg.norm(W - W_last)

                if norm_diff < tolerance:
                    break
                else:
                    # Only update W_last when NOT converged, matching MATLAB logic
                    W_last = W.copy()

                    # TODO: Why does this line exist
            Con[de] = np.array(Con[de])

            Rx = W @ W.conj().T
            CRB = np.trace(np.linalg.inv(FIM))
            # TODO: Why the dot of the first T_k. is omitted
            SR = np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W)))))

            print(f"Configuration: lin={lin[de]}")
            print(f"Converged after {len(Con[de])} iterations")
            print(f"Final objective value: {Con[de][-1, 2]:.6f}")
            print(f"Sum Rate: {SR:.6f}")
            print(f"CRB: {CRB:.6f}")
            print(f"Time elapsed: {perf_counter() - start_time:.6f} seconds")
            print()

            # Plotting
            plt.figure(1)
            plt.plot(Con[de][:, 0], label=f'lin={lin[de]}')
            plt.xlabel('iteration')
            plt.ylabel('Sum Rate (nat/s/Hz)')
            plt.legend(loc='lower right')

            plt.figure(2)
            plt.plot(Con[de][:, 1], label=f'lin={lin[de]}')
            plt.xlabel('iteration')
            plt.ylabel('trace of inverse of FIM')
            plt.legend(loc='lower right')

            plt.figure(3)
            plt.plot(Con[de][:, 2], label=f'lin={lin[de]}')
            plt.xlabel('iteration')
            plt.ylabel('total objective value')
            plt.legend(loc='lower right')

    # Save results - Con should be saved directly to match MATLAB format
    # Convert to object array format that MATLAB uses for cell arrays
    # Create a (3,1) object array to match MATLAB's cell array structure
    Con_matlab = np.empty((3, 1), dtype=object)
    for i in range(3):
        Con_matlab[i, 0] = Con[i]

    # Save to output directory
    output_file = os.path.join(args.output_dir, 'data_convergence.mat')
    os.makedirs(args.output_dir, exist_ok=True)

    savemat(output_file, {
        'Con': Con_matlab  # Use (3,1) object array to match MATLAB cell array
    })
    print(f"Results saved to: {output_file}")

    if args.save_plots:
        # Save plots as files
        plot_names = [
            'fig1_sum_rate_vs_iteration.png',
            'fig1_inverse_fim_trace_vs_iteration.png', 
            'fig1_objective_value_vs_iteration.png'
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

if __name__ == "__main__":
    main()
