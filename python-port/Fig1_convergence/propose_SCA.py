import numpy as np
from scipy import linalg
from scipy.io import savemat
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt
from time import perf_counter
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.steering_matrix import construct_steer_matrix_and_derivative_steer_matrix
from utils.calculate_fim import calculateFIM
from utils.construct_matrixQ import construct_matrixQ

def square_abs(x):
    """Helper function to compute squared absolute value of complex numbers"""
    return np.abs(x) ** 2

def db2pow(x):
    """Convert dB to power"""
    return 10 ** (x / 10)

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
    # Initialize parameters
    Nth = 4
    Ntv = 4
    Nt = Nth * Ntv

    Nrh = 5
    Nrv = 4
    Nr = Nrh * Nrv

    channel_number = 50
    K = 4
    M = 2

    rng = np.random.default_rng(0)
    theta = -np.pi/3 + 2*np.pi/3 * rng.random(M)
    phi = -np.pi/3 + 2*np.pi/3 * rng.random(M)

    num_sensing_streams = Nt
    tolerance = 1e-5

    delta_all = np.arange(-7, 4.8+ 0.01, 0.05)
    Per_all = np.zeros((2, len(delta_all)))
    Con = [[], [], []]

    for channel in [2]:
        for de in range(3):
            lin = [0.05, 0.1, 0.15]
            delta_s = 1
            delta_c = lin[de]

            Pt = db2pow(10-30)  # dBm
            noise_c = db2pow(0-30)  # dBm
            noise_s = db2pow(0-30)  # dBm

            # np.random.seed(channel)
            rng = np.random.default_rng(channel)
            alpha = 0.1 * (1 + 0.2 * rng.standard_normal(M)) * np.exp(1j * 2 * np.pi * rng.random(M))

            L = 30
            kappa = 2 * L / noise_s

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

            for count in range(2000):
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

            print(f"Time elapsed: {perf_counter() - start_time:.6f} seconds")

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
    
    savemat('data_convergence.mat', {
        'Con': Con_matlab  # Use (3,1) object array to match MATLAB cell array
    })
    plt.show()

if __name__ == "__main__":
    main()
