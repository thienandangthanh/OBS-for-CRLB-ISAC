#!/usr/bin/env python3
"""
Steering matrix utilities for UPA (Uniform Planar Array)

This module contains functions for constructing steering matrices and their derivatives
for uniform planar arrays, which are commonly used in MIMO and beamforming applications.
"""

import numpy as np

def construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Mx, My):
    """
    Construct steering matrix and its derivatives for UPA (Uniform Planar Array)

    This function computes the steering matrix for a uniform planar array and its
    derivatives with respect to elevation (theta) and azimuth (phi) angles.

    Args:
        theta (array-like): Elevation angles in radians. Shape: (M,)
        phi (array-like): Azimuth angles in radians. Shape: (M,)
        Mx (int): Number of antenna elements along x-axis
        My (int): Number of antenna elements along y-axis

    Returns:
        tuple: A tuple containing:
            - A (ndarray): Steering matrix of shape (Mx*My, M)
            - dAtheta (ndarray): Derivative of A w.r.t. theta of shape (Mx*My, M)
            - dAphi (ndarray): Derivative of A w.r.t. phi of shape (Mx*My, M)

    Notes:
        - The steering matrix is constructed using the Kronecker product of 
          x-axis and y-axis steering vectors
        - The array elements are assumed to be spaced half-wavelength apart
        - The coordinate system follows the standard spherical coordinate convention

    Example:
        >>> theta = np.array([-np.pi/6, np.pi/6])
        >>> phi = np.array([-np.pi/3, np.pi/4])
        >>> A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(
        ...     theta, phi, 4, 4)
        >>> print(f"Steering matrix shape: {A.shape}")
        Steering matrix shape: (16, 2)
    """
    # Convert inputs to numpy arrays
    theta = np.asarray(theta)
    phi = np.asarray(phi)

    # Generate antenna element indices
    ix = np.arange(Mx).reshape(-1, 1)  # x-axis indices as column vector
    iy = np.arange(My).reshape(-1, 1)  # y-axis indices as column vector

    # Number of steering directions
    M = len(theta)

    # Preallocate output matrices
    A = np.zeros((Mx * My, M), dtype=complex)
    dAtheta = np.zeros((Mx * My, M), dtype=complex)
    dAphi = np.zeros((Mx * My, M), dtype=complex)

    # Compute steering vectors for each direction
    for m in range(M):
        # Steering vectors for x and y dimensions
        # Note: Using pi instead of 2*pi because element spacing is lambda/2
        ax = (1/np.sqrt(Mx)) * np.exp(1j * np.pi * ix * np.sin(theta[m]) * np.sin(phi[m]))
        ay = (1/np.sqrt(My)) * np.exp(1j * np.pi * iy * np.cos(phi[m]))

        # Derivatives of steering vectors
        daxtheta = 1j * np.pi * ix * np.cos(theta[m]) * np.sin(phi[m]) * ax
        daxphi = 1j * np.pi * ix * np.sin(theta[m]) * np.cos(phi[m]) * ax
        dayphi = -1j * np.pi * iy * np.sin(phi[m]) * ay

        # Compute full UPA steering vector using Kronecker product
        A[:, m] = np.kron(ay, ax).flatten()
        dAtheta[:, m] = np.kron(ay, daxtheta).flatten()
        dAphi[:, m] = np.kron(ay, daxphi).flatten() + np.kron(dayphi, ax).flatten()

    return A, dAtheta, dAphi

def validate_steering_matrix_inputs(theta, phi, Mx, My):
    """
    Validate inputs for steering matrix construction

    Args:
        theta (array-like): Elevation angles
        phi (array-like): Azimuth angles  
        Mx (int): Number of elements along x-axis
        My (int): Number of elements along y-axis

    Raises:
        ValueError: If inputs are invalid
    """
    theta = np.asarray(theta)
    phi = np.asarray(phi)

    if theta.shape != phi.shape:
        raise ValueError(f"theta and phi must have the same shape. Got {theta.shape} and {phi.shape}")

    if not isinstance(Mx, int) or Mx <= 0:
        raise ValueError(f"Mx must be a positive integer. Got {Mx}")

    if not isinstance(My, int) or My <= 0:
        raise ValueError(f"My must be a positive integer. Got {My}")

    if np.any(np.abs(theta) > np.pi/2):
        raise ValueError("Elevation angles (theta) should be in range [-π/2, π/2]")

    if np.any(np.abs(phi) > np.pi):
        raise ValueError("Azimuth angles (phi) should be in range [-π, π]")

def construct_steer_matrix_and_derivative_steer_matrix_validated(theta, phi, Mx, My):
    """
    Construct steering matrix with input validation

    This is a wrapper around construct_steer_matrix_and_derivative_steer_matrix
    that includes input validation.

    Args:
        theta (array-like): Elevation angles in radians
        phi (array-like): Azimuth angles in radians
        Mx (int): Number of antenna elements along x-axis
        My (int): Number of antenna elements along y-axis

    Returns:
        tuple: Same as construct_steer_matrix_and_derivative_steer_matrix

    Raises:
        ValueError: If inputs are invalid
    """
    validate_steering_matrix_inputs(theta, phi, Mx, My)
    return construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Mx, My)
