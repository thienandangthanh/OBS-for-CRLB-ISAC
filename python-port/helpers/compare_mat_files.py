#!/usr/bin/env python3
"""
Compare two .mat files for verification of MATLAB to Python translations
Provides detailed comparison reports with numerical tolerance checking
"""

import numpy as np
import scipy.io
import argparse
import sys
import os
from pathlib import Path
import pandas as pd

def load_mat_file(file_path):
    """
    Load .mat file and return data dictionary with metadata removed

    Args:
        file_path (str): Path to .mat file

    Returns:
        dict: Dictionary of variables from .mat file
    """
    try:
        data = scipy.io.loadmat(file_path)
        # Remove MATLAB metadata variables
        clean_data = {k: v for k, v in data.items() if not k.startswith('__')}
        return clean_data
    except Exception as e:
        raise Exception(f"Failed to load {file_path}: {str(e)}")

def compare_arrays(arr1, arr2, var_name, tolerance_abs=1e-12, tolerance_rel=1e-9):
    """
    Compare two numpy arrays with detailed statistics

    Args:
        arr1, arr2: Arrays to compare
        var_name (str): Variable name for reporting
        tolerance_abs (float): Absolute tolerance
        tolerance_rel (float): Relative tolerance

    Returns:
        dict: Comparison results
    """
    result = {
        'variable': var_name,
        'passed': False,
        'shape_match': False,
        'type_match': False,
        'max_abs_error': None,
        'mean_abs_error': None,
        'max_rel_error': None,
        'mean_rel_error': None,
        'details': []
    }

    # Check shapes
    if arr1.shape != arr2.shape:
        result['details'].append(f"Shape mismatch: {arr1.shape} vs {arr2.shape}")
        return result

    result['shape_match'] = True

    # Check data types
    if arr1.dtype != arr2.dtype:
        result['details'].append(f"Data type mismatch: {arr1.dtype} vs {arr2.dtype}")
        # Continue comparison anyway, but note the difference
    else:
        result['type_match'] = True

    # Handle complex numbers
    if np.iscomplexobj(arr1) or np.iscomplexobj(arr2):
        # Ensure both are complex for comparison
        arr1_complex = arr1.astype(complex) if not np.iscomplexobj(arr1) else arr1
        arr2_complex = arr2.astype(complex) if not np.iscomplexobj(arr2) else arr2

        # Compare real and imaginary parts separately
        real_diff = np.abs(arr1_complex.real - arr2_complex.real)
        imag_diff = np.abs(arr1_complex.imag - arr2_complex.imag)
        abs_diff = np.maximum(real_diff, imag_diff)
    else:
        # Real numbers
        abs_diff = np.abs(arr1 - arr2)

    # Calculate statistics
    result['max_abs_error'] = np.max(abs_diff)
    result['mean_abs_error'] = np.mean(abs_diff)

    # Relative error (avoid division by zero)
    arr1_magnitude = np.abs(arr1)
    arr2_magnitude = np.abs(arr2)
    max_magnitude = np.maximum(arr1_magnitude, arr2_magnitude)

    # Only calculate relative error where magnitude is significant
    significant_mask = max_magnitude > tolerance_abs
    if np.any(significant_mask):
        rel_diff = np.zeros_like(abs_diff)
        rel_diff[significant_mask] = abs_diff[significant_mask] / max_magnitude[significant_mask]
        result['max_rel_error'] = np.max(rel_diff)
        result['mean_rel_error'] = np.mean(rel_diff[significant_mask])
    else:
        result['max_rel_error'] = 0.0
        result['mean_rel_error'] = 0.0

    # Check if comparison passes
    abs_pass = result['max_abs_error'] <= tolerance_abs
    rel_pass = result['max_rel_error'] <= tolerance_rel
    result['passed'] = abs_pass or rel_pass

    # Add tolerance information
    result['details'].append(f"Absolute tolerance: {tolerance_abs:.2e}")
    result['details'].append(f"Relative tolerance: {tolerance_rel:.2e}")
    result['details'].append(f"Absolute test: {'PASS' if abs_pass else 'FAIL'}")
    result['details'].append(f"Relative test: {'PASS' if rel_pass else 'FAIL'}")

    return result

def compare_cell_arrays(cell1, cell2, var_name, tolerance_abs=1e-12, tolerance_rel=1e-9):
    """
    Compare MATLAB cell arrays (stored as object arrays in Python)

    Args:
        cell1, cell2: Cell arrays to compare
        var_name (str): Variable name for reporting
        tolerance_abs (float): Absolute tolerance
        tolerance_rel (float): Relative tolerance

    Returns:
        list: List of comparison results for each cell
    """
    results = []

    if cell1.shape != cell2.shape:
        result = {
            'variable': f"{var_name}_shape",
            'passed': False,
            'details': [f"Cell array shape mismatch: {cell1.shape} vs {cell2.shape}"]
        }
        results.append(result)
        return results

    # Compare each cell
    for i, (c1, c2) in enumerate(zip(cell1.flatten(), cell2.flatten())):
        cell_var_name = f"{var_name}_cell_{i}"

        if hasattr(c1, 'shape') and hasattr(c2, 'shape'):
            # Both are arrays
            cell_result = compare_arrays(c1, c2, cell_var_name, tolerance_abs, tolerance_rel)
            results.append(cell_result)
        else:
            # Handle non-array cells
            result = {
                'variable': cell_var_name,
                'passed': np.array_equal(c1, c2),
                'details': [f"Cell content: {c1} vs {c2}"]
            }
            results.append(result)

    return results

def compare_mat_files(file1_path, file2_path, tolerance_abs=1e-12, tolerance_rel=1e-9, 
                      output_report=None, verbose=True):
    """
    Compare two .mat files and generate detailed report

    Args:
        file1_path (str): Path to first .mat file (reference)
        file2_path (str): Path to second .mat file (comparison)
        tolerance_abs (float): Absolute tolerance for numerical comparison
        tolerance_rel (float): Relative tolerance for numerical comparison
        output_report (str): Path to save detailed report (optional)
        verbose (bool): Print detailed output

    Returns:
        dict: Summary of comparison results
    """

    if verbose:
        print(f"Comparing .mat files:")
        print(f"  Reference: {file1_path}")
        print(f"  Comparison: {file2_path}")
        print(f"  Absolute tolerance: {tolerance_abs:.2e}")
        print(f"  Relative tolerance: {tolerance_rel:.2e}")
        print("=" * 80)

    # Load both files
    try:
        data1 = load_mat_file(file1_path)
        data2 = load_mat_file(file2_path)
    except Exception as e:
        print(f"Error loading files: {e}")
        return {'success': False, 'error': str(e)}

    # Get all variable names
    vars1 = set(data1.keys())
    vars2 = set(data2.keys())
    common_vars = vars1.intersection(vars2)
    only_in_file1 = vars1 - vars2
    only_in_file2 = vars2 - vars1

    if verbose:
        print(f"Variables in reference file: {len(vars1)}")
        print(f"Variables in comparison file: {len(vars2)}")
        print(f"Common variables: {len(common_vars)}")
        if only_in_file1:
            print(f"Only in reference: {sorted(only_in_file1)}")
        if only_in_file2:
            print(f"Only in comparison: {sorted(only_in_file2)}")
        print()

    # Compare common variables
    comparison_results = []
    passed_count = 0
    failed_count = 0

    for var_name in sorted(common_vars):
        var1 = data1[var_name]
        var2 = data2[var_name]

        if verbose:
            print(f"Comparing variable: {var_name}")

        # Handle different data types
        if isinstance(var1, np.ndarray) and isinstance(var2, np.ndarray):
            if var1.dtype == np.dtype('O') and var2.dtype == np.dtype('O'):
                # Cell arrays
                cell_results = compare_cell_arrays(var1, var2, var_name, tolerance_abs, tolerance_rel)
                comparison_results.extend(cell_results)

                # Count passes/fails for cells
                for result in cell_results:
                    if result['passed']:
                        passed_count += 1
                    else:
                        failed_count += 1

                if verbose:
                    for result in cell_results:
                        status = "✓ PASS" if result['passed'] else "✗ FAIL"
                        print(f"  {result['variable']}: {status}")
                        if not result['passed']:
                            for detail in result.get('details', []):
                                print(f"    {detail}")
            else:
                # Regular arrays
                result = compare_arrays(var1, var2, var_name, tolerance_abs, tolerance_rel)
                comparison_results.append(result)

                if result['passed']:
                    passed_count += 1
                else:
                    failed_count += 1

                if verbose:
                    status = "✓ PASS" if result['passed'] else "✗ FAIL"
                    print(f"  {status}")
                    if result.get('max_abs_error') is not None:
                        print(f"    Max absolute error: {result['max_abs_error']:.2e}")
                        print(f"    Mean absolute error: {result['mean_abs_error']:.2e}")
                        print(f"    Max relative error: {result['max_rel_error']:.2e}")
                        print(f"    Mean relative error: {result['mean_rel_error']:.2e}")

                    if not result['passed']:
                        for detail in result.get('details', []):
                            print(f"    {detail}")
        else:
            # Non-array comparison
            result = {
                'variable': var_name,
                'passed': var1 == var2,
                'details': [f"Value comparison: {var1} vs {var2}"]
            }
            comparison_results.append(result)

            if result['passed']:
                passed_count += 1
            else:
                failed_count += 1

            if verbose:
                status = "✓ PASS" if result['passed'] else "✗ FAIL"
                print(f"  {status}")
                if not result['passed']:
                    print(f"    {result['details'][0]}")

        if verbose:
            print()

    # Summary
    total_comparisons = passed_count + failed_count
    success_rate = (passed_count / total_comparisons * 100) if total_comparisons > 0 else 0

    summary = {
        'success': True,
        'total_variables': len(common_vars),
        'total_comparisons': total_comparisons,
        'passed': passed_count,
        'failed': failed_count,
        'success_rate': success_rate,
        'only_in_reference': list(only_in_file1),
        'only_in_comparison': list(only_in_file2),
        'detailed_results': comparison_results
    }

    if verbose:
        print("=" * 80)
        print("COMPARISON SUMMARY")
        print("=" * 80)
        print(f"Total variables compared: {len(common_vars)}")
        print(f"Total comparisons performed: {total_comparisons}")
        print(f"Passed: {passed_count}")
        print(f"Failed: {failed_count}")
        print(f"Success rate: {success_rate:.1f}%")

        if failed_count == 0:
            print("\n🎉 All comparisons PASSED! Files are equivalent within tolerance.")
        else:
            print(f"\n⚠️  {failed_count} comparison(s) FAILED. Check details above.")

    # Save detailed report if requested
    if output_report:
        save_comparison_report(summary, file1_path, file2_path, output_report)
        if verbose:
            print(f"\nDetailed report saved to: {output_report}")

    return summary

def save_comparison_report(summary, file1_path, file2_path, output_path):
    """
    Save detailed comparison report to Excel file

    Args:
        summary (dict): Comparison summary
        file1_path (str): Path to first file
        file2_path (str): Path to second file
        output_path (str): Path to save report
    """

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        # Summary sheet
        summary_data = {
            'Metric': ['Reference File', 'Comparison File', 'Total Variables', 
                       'Total Comparisons', 'Passed', 'Failed', 'Success Rate (%)'],
            'Value': [file1_path, file2_path, summary['total_variables'],
                      summary['total_comparisons'], summary['passed'], 
                      summary['failed'], f"{summary['success_rate']:.1f}"]
        }

        df_summary = pd.DataFrame(summary_data)
        df_summary.to_excel(writer, sheet_name='Summary', index=False)

        # Detailed results sheet
        if summary['detailed_results']:
            results_data = []
            for result in summary['detailed_results']:
                row = {
                    'Variable': result['variable'],
                    'Status': 'PASS' if result['passed'] else 'FAIL',
                    'Shape Match': result.get('shape_match', 'N/A'),
                    'Type Match': result.get('type_match', 'N/A'),
                    'Max Abs Error': result.get('max_abs_error', 'N/A'),
                    'Mean Abs Error': result.get('mean_abs_error', 'N/A'),
                    'Max Rel Error': result.get('max_rel_error', 'N/A'),
                    'Mean Rel Error': result.get('mean_rel_error', 'N/A'),
                    'Details': '; '.join(result.get('details', []))
                }
                results_data.append(row)

            df_results = pd.DataFrame(results_data)
            df_results.to_excel(writer, sheet_name='Detailed_Results', index=False)

        # Variables only in one file
        if summary['only_in_reference'] or summary['only_in_comparison']:
            unique_vars_data = []

            for var in summary['only_in_reference']:
                unique_vars_data.append({'Variable': var, 'Location': 'Reference Only'})

            for var in summary['only_in_comparison']:
                unique_vars_data.append({'Variable': var, 'Location': 'Comparison Only'})

            if unique_vars_data:
                df_unique = pd.DataFrame(unique_vars_data)
                df_unique.to_excel(writer, sheet_name='Unique_Variables', index=False)

def main():
    parser = argparse.ArgumentParser(
        description='Compare two .mat files for verification of MATLAB to Python translations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic comparison
  python compare_mat_files.py matlab_results.mat python_results.mat

  # With custom tolerances
  python compare_mat_files.py file1.mat file2.mat --abs-tol 1e-10 --rel-tol 1e-8

  # Save detailed report
  python compare_mat_files.py file1.mat file2.mat --report comparison_report.xlsx

  # Quiet mode (only show summary)
  python compare_mat_files.py file1.mat file2.mat --quiet
        """
    )

    parser.add_argument('file1', help='Path to first .mat file (reference)')
    parser.add_argument('file2', help='Path to second .mat file (comparison)')

    parser.add_argument('--abs-tol', type=float, default=1e-12,
                        help='Absolute tolerance for numerical comparison (default: 1e-12)')
    parser.add_argument('--rel-tol', type=float, default=1e-9,
                        help='Relative tolerance for numerical comparison (default: 1e-9)')
    parser.add_argument('--report', type=str,
                        help='Save detailed comparison report to Excel file')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress detailed output, show only summary')

    args = parser.parse_args()

    # Validate input files
    for file_path in [args.file1, args.file2]:
        if not os.path.exists(file_path):
            print(f"Error: File '{file_path}' does not exist!")
            sys.exit(1)
        if not file_path.lower().endswith('.mat'):
            print(f"Warning: File '{file_path}' does not have .mat extension")

    # Perform comparison
    try:
        summary = compare_mat_files(
            args.file1, 
            args.file2,
            tolerance_abs=args.abs_tol,
            tolerance_rel=args.rel_tol,
            output_report=args.report,
            verbose=not args.quiet
        )

        if not summary['success']:
            sys.exit(1)

        # Exit with error code if comparisons failed
        if summary['failed'] > 0:
            sys.exit(1)

    except Exception as e:
        print(f"Error during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
