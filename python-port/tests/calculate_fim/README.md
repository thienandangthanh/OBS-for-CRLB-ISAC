# Fisher Information Matrix (FIM) Tests - MATLAB Reference Validation

This directory contains comprehensive tests for validating the Python implementation of the `calculateFIM` function against the MATLAB reference implementation.

## Overview

The Fisher Information Matrix (FIM) is a crucial component in Integrated Sensing and Communication (ISAC) systems for:
- Parameter estimation accuracy analysis
- Cramer-Rao Lower Bound (CRLB) computation
- Beamforming optimization

## Test Approach

### **MATLAB as Reference Implementation**
- MATLAB `calculateFIM` is the **trusted reference implementation**
- Generates canonical test data with known-good results
- Saves artifacts in a location-independent manner

### **Python as Implementation Under Test**
- Python implementation validates against MATLAB reference data
- Focuses on proving equivalence with the reference
- Reports differences and validation status

## Test Structure

### MATLAB Reference Generation (`OBS_for_CRLB_ISAC/tests/`)
- **`test_calculate_fim_generate_data.m`**: Generates reference test data using MATLAB `calculateFIM`
- **`test_calculate_fim_data.mat`**: Generated reference data file (saved next to script)

### Python Validation (`python-port/tests/calculate_fim/`)
- **`validate_against_matlab.py`**: Validates Python implementation against MATLAB reference

## Running the Tests

### **Step 1: Generate MATLAB Reference Data**
```matlab
% From any directory in MATLAB
test_calculate_fim_generate_data()
```

This will:
1. Generate test parameters using fixed random seed (reproducible)
2. Compute FIM using MATLAB `calculateFIM` function
3. Validate MATLAB results for consistency
4. Save reference data to `OBS_for_CRLB_ISAC/tests/test_calculate_fim_data.mat`

### **Step 2: Validate Python Implementation**
```bash
# From project root or calculate_fim test directory
cd python-port/tests/calculate_fim
python validate_against_matlab.py
```

This will:
1. Load MATLAB reference data automatically
2. Compute FIM using Python `calculate_fim` and `calculateFIM`
3. Compare Python results vs MATLAB reference
4. Report validation status and differences

## Test Parameters

- **Array Configuration**: 4×4 transmit, 5×4 receive antenna arrays
- **Targets**: 2 targets with random angles and reflection coefficients
- **Beamforming**: Communication + sensing streams
- **Noise**: Configurable sensing noise level (1e-3 default)
- **Snapshots**: L=30 symbols
- **Random Seed**: Fixed (42) for reproducibility

## Validation Criteria

### Numerical Accuracy
- **Absolute tolerance**: 1e-10 (relaxed for numerical precision)
- **Relative tolerance**: 1e-8 (accounts for small denominators)
- **Machine precision override**: Auto-pass for differences < 1e-12

### Mathematical Properties
- **Hermitian**: FIM = FIM†
- **Positive semi-definite**: All eigenvalues ≥ 0
- **Scaling**: Linear with L, inverse with noise_s
- **Input equivalence**: W matrix vs Rx covariance matrix

## Output Files

**`test_calculate_fim_data.mat`**: MATLAB-generated reference data

## Expected Results

When properly implemented, validation should show:
- ✅ Maximum differences at machine precision level (< 1e-12)
- ✅ Both implementations produce Hermitian, positive semi-definite matrices
- ✅ Identical mathematical properties
- ✅ Correct scaling behavior

## Key Functions Tested

### Python (`utils/calculate_fim.py`)
- `calculate_fim()`: Main function supporting W or Rx inputs
- `calculateFIM()`: Legacy wrapper for backward compatibility

### MATLAB (`utils/calculateFIM.m`)
- `calculateFIM()`: Reference implementation

## Troubleshooting

### Common Issues
1. **MATLAB data not found**: Run `test_calculate_fim_generate_data()` first
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