# Matrix Q Construction Tests - MATLAB Reference Validation

This directory contains comprehensive tests for validating the Python implementation of the `construct_matrixQ` function against the MATLAB reference implementation.

## Overview

The matrix Q construction is a crucial component in Integrated Sensing and Communication (ISAC) systems for:
- Beamforming optimization
- Fisher Information Matrix (FIM) computation
- Cramer-Rao Lower Bound (CRLB) analysis

## Test Approach

### **MATLAB as Reference Implementation**
- MATLAB `construct_matrixQ` is the **trusted reference implementation**
- Generates canonical test data with known-good results
- Saves artifacts in a location-independent manner

### **Python as Implementation Under Test**
- Python implementation validates against MATLAB reference data
- Focuses on proving equivalence with the reference
- Reports differences and validation status

## Test Structure

### MATLAB Reference Generation (`OBS_for_CRLB_ISAC/tests/`)
- **`test_construct_matrixQ_generate_data.m`**: Generates reference test data using MATLAB `construct_matrixQ`
- **`test_construct_matrixQ_data.mat`**: Generated reference data file (saved next to script)

### Python Validation (`python-port/tests/construct_matrixQ/`)
- **`validate_against_matlab.py`**: Validates Python implementation against MATLAB reference

## Running the Tests

### **Step 1: Generate MATLAB Reference Data**
```matlab
% From any directory in MATLAB
test_construct_matrixQ_generate_data()
```

This will:
1. Generate test parameters using fixed random seed (reproducible)
2. Compute Q matrix using MATLAB `construct_matrixQ` function
3. Validate MATLAB results for consistency
4. Save reference data to `OBS_for_CRLB_ISAC/tests/test_construct_matrixQ_data.mat`

### **Step 2: Validate Python Implementation**
```bash
# From project root or construct_matrixQ test directory
cd python-port/tests/construct_matrixQ
python validate_against_matlab.py
```

This will:
1. Load MATLAB reference data automatically
2. Compute Q matrix using Python `construct_matrixQ`
3. Compare Python results vs MATLAB reference
4. Report validation status and differences

## Test Parameters

The test suite includes three test cases with varying parameters:
1. **Small matrices**: Basic test case with small dimensions
2. **Larger matrices**: Test case with larger dimensions and different parameters
3. **Edge case**: Test with very small values to check numerical stability

## Validation Criteria

### Numerical Accuracy
- **Absolute tolerance**: 1e-6 (realistic for multiple matrix multiplications)
- **Relative tolerance**: 1e-8 (accounts for small denominators)
- **Machine precision override**: Auto-pass for differences < 1e-12

### Mathematical Properties
- **Square matrix**: Q must be square
- **Scaling**: Linear with L, inverse with noise_s

## Expected Results

When properly implemented, validation should show:
- ✅ Maximum differences within numerical tolerance (< 1e-6 absolute, < 1e-8 relative)
- ✅ Both implementations produce matrices with identical properties
- ✅ Correct scaling behavior

Note: The Q matrix is not guaranteed to be Hermitian in the general case, matching MATLAB's behavior.

## Key Functions Tested

### Python (`utils/construct_matrixQ.py`)
- `construct_matrixQ()`: Main function for Q matrix construction

### MATLAB (`utils/construct_matrixQ.m`)
- `construct_matrixQ()`: Reference implementation

## Troubleshooting

### Common Issues
1. **MATLAB data not found**: Run `test_construct_matrixQ_generate_data()` first
2. **Path errors**: Scripts auto-detect relative paths
3. **Data loading errors**: Check MATLAB data structure format

### Debug Tips
- Use `validate_against_matlab.py` for detailed comparison
- Check intermediate matrix properties
- Verify random seed consistency between MATLAB/Python
- Compare computation times as sanity check

## Implementation Notes

### MATLAB Reference Generation
- Uses `mfilename('fullpath')` for location-independent saving
- Fixed random seed for reproducibility
- Comprehensive validation of reference data quality

### Python Validation
- Robust data loading with nested structure handling
- Both absolute and relative difference analysis
- Automatic tolerance assessment with overrides
- Mathematical property verification

This approach provides a robust, hierarchical validation framework where MATLAB serves as the authoritative reference implementation. 