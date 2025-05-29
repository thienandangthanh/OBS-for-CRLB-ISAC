"""
Tests for construct_matrixQ function.

This package contains tests to validate the Python implementation of the
matrix Q construction against the MATLAB reference implementation.
"""

from .validate_against_matlab import MATLABReferenceValidator, main

__all__ = ['MATLABReferenceValidator', 'main'] 