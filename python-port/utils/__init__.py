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

from .calculate_fim import (
    calculate_fim,
    calculateFIM
)

from .construct_matrixQ import (
    construct_matrixQ
)

from .simulation_config import (
    SimulationConfig,
    db2pow
)

__all__ = [
    'construct_steer_matrix_and_derivative_steer_matrix',
    'construct_steer_matrix_and_derivative_steer_matrix_validated', 
    'validate_steering_matrix_inputs',
    'calculate_fim',
    'calculateFIM',
    'construct_matrixQ',
    'SimulationConfig',
    'db2pow'
]
