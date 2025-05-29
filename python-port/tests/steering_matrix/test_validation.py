#!/usr/bin/env python3
"""
Test script for steering matrix input validation
"""

import numpy as np
import pytest
import sys
import os

# Add the parent directories to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from utils.steering_matrix import (
    construct_steer_matrix_and_derivative_steer_matrix_validated,
    validate_steering_matrix_inputs
)

def test_valid_inputs():
    """Test with valid inputs"""
    theta = np.array([-np.pi/6, np.pi/6])
    phi = np.array([-np.pi/3, np.pi/4])
    Mx, My = 4, 4
    
    # Should not raise any exception
    A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix_validated(
        theta, phi, Mx, My
    )
    
    assert A.shape == (Mx * My, len(theta))
    assert dAtheta.shape == (Mx * My, len(theta))
    assert dAphi.shape == (Mx * My, len(theta))

def test_mismatched_angle_shapes():
    """Test with mismatched theta and phi shapes"""
    theta = np.array([-np.pi/6, np.pi/6])
    phi = np.array([-np.pi/3])  # Different length
    Mx, My = 4, 4
    
    with pytest.raises(ValueError, match="theta and phi must have the same shape"):
        validate_steering_matrix_inputs(theta, phi, Mx, My)

def test_invalid_array_dimensions():
    """Test with invalid array dimensions"""
    theta = np.array([-np.pi/6, np.pi/6])
    phi = np.array([-np.pi/3, np.pi/4])
    
    # Test negative Mx
    with pytest.raises(ValueError, match="Mx must be a positive integer"):
        validate_steering_matrix_inputs(theta, phi, -1, 4)
    
    # Test zero My
    with pytest.raises(ValueError, match="My must be a positive integer"):
        validate_steering_matrix_inputs(theta, phi, 4, 0)
    
    # Test non-integer Mx
    with pytest.raises(ValueError, match="Mx must be a positive integer"):
        validate_steering_matrix_inputs(theta, phi, 4.5, 4)

def test_angle_range_validation():
    """Test angle range validation"""
    Mx, My = 4, 4
    
    # Test theta out of range
    theta = np.array([np.pi])  # > π/2
    phi = np.array([0])
    with pytest.raises(ValueError, match="Elevation angles.*should be in range"):
        validate_steering_matrix_inputs(theta, phi, Mx, My)
    
    # Test phi out of range
    theta = np.array([0])
    phi = np.array([2*np.pi])  # > π
    with pytest.raises(ValueError, match="Azimuth angles.*should be in range"):
        validate_steering_matrix_inputs(theta, phi, Mx, My)

def test_edge_cases():
    """Test edge cases"""
    # Single angle
    theta = np.array([0])
    phi = np.array([0])
    Mx, My = 2, 2
    
    A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix_validated(
        theta, phi, Mx, My
    )
    
    assert A.shape == (4, 1)
    assert dAtheta.shape == (4, 1)
    assert dAphi.shape == (4, 1)
    
    # Minimum array size
    Mx, My = 1, 1
    A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix_validated(
        theta, phi, Mx, My
    )
    
    assert A.shape == (1, 1)
    assert dAtheta.shape == (1, 1)
    assert dAphi.shape == (1, 1)

if __name__ == "__main__":
    # Run tests manually if pytest is not available
    try:
        test_valid_inputs()
        test_mismatched_angle_shapes()
        test_invalid_array_dimensions()
        test_angle_range_validation()
        test_edge_cases()
        print("All validation tests passed!")
    except Exception as e:
        print(f"Test failed: {e}")
        sys.exit(1)
