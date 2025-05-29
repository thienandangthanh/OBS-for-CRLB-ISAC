"""
Construct matrix Q for optimization

This module provides the construct_matrixQ function that constructs
matrix Q used in the optimization process for ISAC systems.
"""

import numpy as np


def construct_matrixQ(L, noise_s, Phi, A, dAtheta, dAphi, B, dBtheta, dBphi, U):
    """
    Construct matrix Q for optimization
    
    Parameters:
    -----------
    L : float
        Number of symbols
    noise_s : float
        Sensing noise power
    Phi : ndarray
        Covariance matrix (CRBM*CRBM)
    A : ndarray
        Transmit steering matrix
    dAtheta : ndarray
        Derivative of A w.r.t. theta
    dAphi : ndarray
        Derivative of A w.r.t. phi
    B : ndarray
        Receive steering matrix
    dBtheta : ndarray
        Derivative of B w.r.t. theta
    dBphi : ndarray
        Derivative of B w.r.t. phi
    U : ndarray
        Diagonal matrix of channel coefficients
    
    Returns:
    --------
    Q : ndarray
        Constructed matrix Q
    """
    M = U.shape[0]

    # Extract Phi blocks
    phi11 = Phi[:M, :M]
    phi12 = Phi[:M, M:2*M]
    phi13 = Phi[:M, 2*M:3*M]
    phi14 = Phi[:M, 3*M:4*M]
    phi22 = Phi[M:2*M, M:2*M]
    phi23 = Phi[M:2*M, 2*M:3*M]
    phi24 = Phi[M:2*M, 3*M:4*M]
    phi33 = Phi[2*M:3*M, 2*M:3*M]
    phi34 = Phi[2*M:3*M, 3*M:4*M]
    phi44 = Phi[3*M:4*M, 3*M:4*M]

    # Compute Q components following MATLAB structure exactly
    # Q11 = (A * U' * (phi11 .* (dBtheta' * dBtheta)) * U * A')  + ...
    Q11 = (A @ U.conj().T @ (phi11 * (dBtheta.conj().T @ dBtheta)) @ U @ A.conj().T) + \
          (dAtheta @ U.conj().T @ (phi11 * (B.conj().T @ dBtheta)) @ U @ A.conj().T) + \
          (A @ U.conj().T @ (phi11 * (dBtheta.conj().T @ B)) @ U @ dAtheta.conj().T) + \
          (dAtheta @ U.conj().T @ (phi11 * (B.conj().T @ B)) @ U @ dAtheta.conj().T)

    # Q12 = (A * U' * (2*phi12.* (dBtheta' * dBphi)) * U * A' )  + ...
    Q12 = (A @ U.conj().T @ (2*phi12 * (dBtheta.conj().T @ dBphi)) @ U @ A.conj().T) + \
          (dAtheta @ U.conj().T @ (2*phi12 * (B.conj().T @ dBphi)) @ U @ A.conj().T) + \
          (A @ U.conj().T @ (2*phi12 * (dBtheta.conj().T @ B)) @ U @ dAphi.conj().T) + \
          (dAtheta @ U.conj().T @ (2*phi12 * (B.conj().T @ B)) @ U @ dAphi.conj().T)

    # Q13 = (A * U' * ((2*phi13+2j*phi14).* (dBtheta' * B)) * A' )  + ...
    Q13 = (A @ U.conj().T @ ((2*phi13+2j*phi14) * (dBtheta.conj().T @ B)) @ A.conj().T) + \
          (dAtheta @ U.conj().T @ ((2*phi13+2j*phi14) * (B.conj().T @ B)) @ A.conj().T)

    # Q22 = (A * U' * (phi22.* (dBphi' * dBphi)) * U * A' )  + ...
    Q22 = (A @ U.conj().T @ (phi22 * (dBphi.conj().T @ dBphi)) @ U @ A.conj().T) + \
          (dAphi @ U.conj().T @ (phi22 * (B.conj().T @ dBphi)) @ U @ A.conj().T) + \
          (A @ U.conj().T @ (phi22 * (dBphi.conj().T @ B)) @ U @ dAphi.conj().T) + \
          (dAphi @ U.conj().T @ (phi22 * (B.conj().T @ B)) @ U @ dAphi.conj().T)

    # Q23 = (A * U' * ((2*phi23+2j*phi24).* (dBphi' * B)) *A')  + ...
    Q23 = (A @ U.conj().T @ ((2*phi23+2j*phi24) * (dBphi.conj().T @ B)) @ A.conj().T) + \
          (dAphi @ U.conj().T @ ((2*phi23+2j*phi24) * (B.conj().T @ B)) @ A.conj().T)

    # Q33 = ( A * ((phi33+phi44+2j*phi34).* (B' * B))*A') ;
    Q33 = A @ ((phi33+phi44+2j*phi34) * (B.conj().T @ B)) @ A.conj().T

    # Q=2*L/noise_s*(Q11+Q12+Q13+Q22+Q23+Q33);
    Q = 2*L/noise_s*(Q11+Q12+Q13+Q22+Q23+Q33)
    
    return Q 