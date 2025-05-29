"""
Fisher Information Matrix (FIM) calculation for ISAC systems.

This module contains functions for computing the Fisher Information Matrix
used in Integrated Sensing and Communication (ISAC) systems for parameter estimation
and Cramer-Rao Lower Bound (CRLB) analysis.
"""

import numpy as np


def calculate_fim(L, noise_s, W_or_Rx, A, dAtheta, dAphi, B, dBtheta, dBphi, U):
    """
    Calculate Fisher Information Matrix for ISAC system.

    This function computes the Fisher Information Matrix (FIM) for 
    Integrated Sensing and Communication (ISAC) systems, used to derive
    the Cramer-Rao Lower Bound (CRLB) for parameter estimation.

    Args:
        L (float): Number of symbols/snapshots
        noise_s (float): Sensing noise power  
        W_or_Rx (np.ndarray): Either beamforming matrix W (Nt x total_streams) 
                             or receive covariance matrix Rx (Nt x Nt)
        A (np.ndarray): Transmit steering matrix (Nt x M)
        dAtheta (np.ndarray): Partial derivative of A w.r.t. elevation angle theta (Nt x M)
        dAphi (np.ndarray): Partial derivative of A w.r.t. azimuth angle phi (Nt x M)
        B (np.ndarray): Receive steering matrix (Nr x M) 
        dBtheta (np.ndarray): Partial derivative of B w.r.t. elevation angle theta (Nr x M)
        dBphi (np.ndarray): Partial derivative of B w.r.t. azimuth angle phi (Nr x M)
        U (np.ndarray): Diagonal matrix of reflection coefficients (M x M)

    Returns:
        np.ndarray: Fisher Information Matrix (4M x 4M)
                   Block structure: [theta_params; phi_params; real_alpha; imag_alpha]

    Note:
        The FIM is constructed for parameter vector:
        [theta_1, ..., theta_M, phi_1, ..., phi_M, Re(alpha_1), ..., Re(alpha_M), Im(alpha_1), ..., Im(alpha_M)]
    """
    # Determine if input is W (beamforming matrix) or Rx (covariance matrix)
    Nt, num_cols = W_or_Rx.shape
    if Nt == num_cols:
        # Square matrix - assume it's Rx (covariance matrix)
        Rx = W_or_Rx
    else:
        # Non-square matrix - assume it's W (beamforming matrix)
        Rx = W_or_Rx @ W_or_Rx.conj().T

    # Get number of targets
    M = U.shape[0]

    # Compute intermediate Fisher information terms
    # F11: theta-theta block
    F11 = (U @ A.conj().T @ Rx @ A @ U.conj().T).T * (dBtheta.conj().T @ dBtheta) + \
        (U @ A.conj().T @ Rx @ dAtheta @ U.conj().T).T * (B.conj().T @ dBtheta) + \
        (U @ dAtheta.conj().T @ Rx @ A @ U.conj().T).T * (dBtheta.conj().T @ B) + \
        (U @ dAtheta.conj().T @ Rx @ dAtheta @ U.conj().T).T * (B.conj().T @ B)

    # F12: theta-phi block  
    F12 = (U @ A.conj().T @ Rx @ A @ U.conj().T).T * (dBtheta.conj().T @ dBphi) + \
        (U @ A.conj().T @ Rx @ dAtheta @ U.conj().T).T * (B.conj().T @ dBphi) + \
        (U @ dAphi.conj().T @ Rx @ A @ U.conj().T).T * (dBtheta.conj().T @ B) + \
        (U @ dAphi.conj().T @ Rx @ dAtheta @ U.conj().T).T * (B.conj().T @ B)

    # F13: theta-alpha block (complex coefficients)
    F13 = (A.conj().T @ Rx @ A @ U.conj().T).T * (dBtheta.conj().T @ B) + \
        (A.conj().T @ Rx @ dAtheta @ U.conj().T).T * (B.conj().T @ B)

    # F22: phi-phi block
    F22 = (U @ A.conj().T @ Rx @ A @ U.conj().T).T * (dBphi.conj().T @ dBphi) + \
        (U @ A.conj().T @ Rx @ dAphi @ U.conj().T).T * (B.conj().T @ dBphi) + \
        (U @ dAphi.conj().T @ Rx @ A @ U.conj().T).T * (dBphi.conj().T @ B) + \
        (U @ dAphi.conj().T @ Rx @ dAphi @ U.conj().T).T * (B.conj().T @ B)

    # F23: phi-alpha block (complex coefficients)  
    F23 = (A.conj().T @ Rx @ A @ U.conj().T).T * (dBphi.conj().T @ B) + \
        (A.conj().T @ Rx @ dAphi @ U.conj().T).T * (B.conj().T @ B)

    # F33: alpha-alpha block (complex coefficients)
    F33 = (A.conj().T @ Rx @ A).T * (B.conj().T @ B)

    # Initialize FIM
    FIM = np.zeros((4*M, 4*M), dtype=complex)

    # Fill the FIM blocks
    # Real parts for angle parameters and real parts of complex coefficients
    FIM[:M, :M] = np.real(F11)                    # theta-theta  
    FIM[:M, M:2*M] = np.real(F12)                 # theta-phi
    FIM[:M, 2*M:3*M] = np.real(F13)               # theta-Re(alpha)
    FIM[:M, 3*M:4*M] = -np.imag(F13)              # theta-Im(alpha)
    FIM[M:2*M, M:2*M] = np.real(F22)              # phi-phi
    FIM[M:2*M, 2*M:3*M] = np.real(F23)            # phi-Re(alpha)  
    FIM[M:2*M, 3*M:4*M] = -np.imag(F23)           # phi-Im(alpha)
    FIM[2*M:3*M, 2*M:3*M] = np.real(F33)          # Re(alpha)-Re(alpha)
    FIM[2*M:3*M, 3*M:4*M] = -np.imag(F33)         # Re(alpha)-Im(alpha)
    FIM[3*M:4*M, 3*M:4*M] = np.real(F33)          # Im(alpha)-Im(alpha)

    # Make the matrix Hermitian by symmetry
    FIM = np.triu(FIM) + np.triu(FIM).conj().T - np.diag(np.diag(FIM))

    # Apply scaling factor
    return 2 * L / noise_s * FIM


# Legacy function name for backward compatibility
def calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U):
    """
    Legacy wrapper for calculate_fim with W parameter only (backward compatibility).

    Args:
        W (np.ndarray): Beamforming matrix (Nt x total_streams)
        (other parameters same as calculate_fim)

    Returns:
        np.ndarray: Fisher Information Matrix (4M x 4M)
    """
    return calculate_fim(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
