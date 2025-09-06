"""
Pytest wrapper for construct_matrixQ validation tests
"""

import pytest
from .validate_against_matlab import MATLABReferenceValidator

@pytest.fixture(scope="module")
def validator():
    """Create validator instance to be used by all tests"""
    return MATLABReferenceValidator()

def test_matlab_reference_validation(validator):
    """Test Python implementation against MATLAB reference data"""
    all_passed, results = validator.validate_all_test_cases()
    assert all_passed, "Validation against MATLAB reference failed"
    
    # Check results for each test case
    for test_case, result in results.items():
        assert result['max_abs_diff'] < 1e-6, f"{test_case}: Maximum absolute difference too high: {result['max_abs_diff']}"
        assert result['max_rel_diff'] < 1e-8, f"{test_case}: Maximum relative difference too high: {result['max_rel_diff']}"

@pytest.mark.parametrize("test_case", ["test1", "test2", "test3"])
def test_mathematical_properties(validator, test_case):
    """Test mathematical properties for each test case"""
    validator.test_mathematical_properties(test_case) 