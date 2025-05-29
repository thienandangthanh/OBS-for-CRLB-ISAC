"""
Pytest wrapper for calculate_fim validation tests
"""

import pytest
from .validate_against_matlab import MATLABReferenceValidator

@pytest.fixture(scope="module")
def validator():
    """Create validator instance to be used by all tests"""
    return MATLABReferenceValidator()

def test_matlab_reference_validation(validator):
    """Test Python implementation against MATLAB reference data"""
    validation_passed, results = validator.validate_python_implementation()
    
    # Use the same validation logic as in validate_against_matlab.py
    max_abs_diff = results['max_abs_diff']
    
    # Test passes if either:
    # 1. Differences are at machine precision level (< 1e-12)
    # 2. Differences are within acceptable tolerance (< 1e-10)
    is_machine_precision = max_abs_diff < 1e-12
    is_within_tolerance = max_abs_diff < 1e-10
    
    assert is_machine_precision or is_within_tolerance, (
        f"Validation failed: max_abs_diff={max_abs_diff:.2e} is neither at machine precision "
        f"(< 1e-12) nor within tolerance (< 1e-10)"
    )

def test_input_variations(validator):
    """Test various input scenarios"""
    validator.test_input_variations()
    
    # Add assertions for input variations
    # These could check specific properties that should hold true
    # For example, scaling properties, matrix properties, etc. 