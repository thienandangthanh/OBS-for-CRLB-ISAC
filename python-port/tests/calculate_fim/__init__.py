"""
Tests for calculateFIM function.

This package contains tests to validate the Python implementation of the
Fisher Information Matrix calculation against the MATLAB reference implementation.
"""

from .validate_against_matlab import MATLABReferenceValidator, main

__all__ = ['MATLABReferenceValidator', 'main'] 