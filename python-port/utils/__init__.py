"""
Utilities package for OBS-for-CRLB-ISAC

This package contains utility functions and modules used across different
parts of the OBS-for-CRLB-ISAC project.
"""

from .steering_matrix import (
    construct_steer_matrix_and_derivative_steer_matrix,
    construct_steer_matrix_and_derivative_steer_matrix_validated,
    validate_steering_matrix_inputs
)

__all__ = [
    'construct_steer_matrix_and_derivative_steer_matrix',
    'construct_steer_matrix_and_derivative_steer_matrix_validated', 
    'validate_steering_matrix_inputs'
]
