# Tests Directory

This directory contains test files and scripts for validating the functionality of the python-port tools.

## Files

### `mod-mat-metadata.py`
- **Purpose**: Generates modified .mat files with altered metadata for testing
- **Usage**: Creates test artifacts to verify that `compare_mat_files.py` correctly handles metadata differences
- **Arguments**: 
  - `input`: Path to input .mat file
  - `output`: Path to output .mat file
- **Example**: `python mod-mat-metadata.py ../Fig1_convergence/data_convergence.mat mod_data_convergence.mat`

### `test_compare_mat_files.py`
- **Purpose**: Automated test for the `compare_mat_files.py` script
- **Usage**: Verifies that metadata differences don't affect data comparison results
- **Dependencies**: Uses both `mod-mat-metadata.py` and `../helpers/compare_mat_files.py`

## Running Tests

### Individual Test
```bash
cd python-port/tests
python test_compare_mat_files.py
```

### Generate Test Data Only
```bash
cd python-port/tests
python mod-mat-metadata.py ../Fig1_convergence/data_convergence.mat mod_data_convergence.mat
```

### Manual Comparison
```bash
cd python-port/tests
python ../helpers/compare_mat_files.py ../Fig1_convergence/data_convergence.mat mod_data_convergence.mat
```

## Expected Results

When running the tests, you should see:
1. **Metadata differences detected**: The files have different metadata (creation time, software version, etc.)
2. **Data comparisons pass**: All actual data variables are identical within tolerance
3. **Test passes**: Confirms that `compare_mat_files.py` correctly separates metadata from data

This validates that the comparison tool is working correctly and won't falsely report data differences when only metadata differs between MATLAB and Python generated files.

## Cleanup

Test files are automatically cleaned up after running `test_compare_mat_files.py`. To manually clean up:
```bash
cd python-port/tests
rm -f mod_data_convergence.mat
```
