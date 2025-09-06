#!/usr/bin/env python3
"""
Test script for construct_steer_matrix_and_derivative_steer_matrix function
Compares Python implementation with MATLAB reference data
"""

import numpy as np
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt
import sys
import os

# Add the parent directories to the path to find the function
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from utils.steering_matrix import construct_steer_matrix_and_derivative_steer_matrix

def compare_complex_matrices(mat1, mat2, tolerance=1e-12):
    """Compare two complex matrices with a tolerance"""
    if mat1.shape != mat2.shape:
        return False, np.inf, "Shape mismatch"
    
    diff = np.abs(mat1 - mat2)
    max_diff = float(np.max(diff))
    mean_diff = float(np.mean(diff))
    is_close = bool(max_diff < tolerance)
    
    return is_close, max_diff, mean_diff

def extract_matlab_data(case_data):
    """Extract data from MATLAB structure, handling different possible formats"""
    try:
        # Try direct access first
        inputs = case_data['inputs'][0, 0]
        outputs = case_data['outputs'][0, 0]
        
        # Extract inputs
        theta = inputs['theta'][0, 0].flatten()
        phi = inputs['phi'][0, 0].flatten()
        Mx = int(inputs['Mx'][0, 0])
        My = int(inputs['My'][0, 0])
        
        # Extract outputs
        A_matlab = outputs['A'][0, 0]
        dAtheta_matlab = outputs['dAtheta'][0, 0]
        dAphi_matlab = outputs['dAphi'][0, 0]
        
        return {
            'inputs': {'theta': theta, 'phi': phi, 'Mx': Mx, 'My': My},
            'outputs': {'A': A_matlab, 'dAtheta': dAtheta_matlab, 'dAphi': dAphi_matlab}
        }
        
    except (KeyError, IndexError, TypeError) as e:
        print(f"Error extracting data (method 1): {e}")
        
        # Try alternative access method
        try:
            # Sometimes MATLAB structures are accessed differently
            inputs = case_data[0, 0]['inputs'][0, 0]
            outputs = case_data[0, 0]['outputs'][0, 0]
            
            theta = inputs['theta'][0, 0].flatten()
            phi = inputs['phi'][0, 0].flatten()
            Mx = int(inputs['Mx'][0, 0])
            My = int(inputs['My'][0, 0])
            
            A_matlab = outputs['A'][0, 0]
            dAtheta_matlab = outputs['dAtheta'][0, 0]
            dAphi_matlab = outputs['dAphi'][0, 0]
            
            return {
                'inputs': {'theta': theta, 'phi': phi, 'Mx': Mx, 'My': My},
                'outputs': {'A': A_matlab, 'dAtheta': dAtheta_matlab, 'dAphi': dAphi_matlab}
            }
            
        except (KeyError, IndexError, TypeError) as e2:
            print(f"Error extracting data (method 2): {e2}")
            
            # Try yet another method - direct field access
            try:
                theta = case_data['inputs']['theta'].flatten()
                phi = case_data['inputs']['phi'].flatten()
                Mx = int(case_data['inputs']['Mx'])
                My = int(case_data['inputs']['My'])
                
                A_matlab = case_data['outputs']['A']
                dAtheta_matlab = case_data['outputs']['dAtheta']
                dAphi_matlab = case_data['outputs']['dAphi']
                
                return {
                    'inputs': {'theta': theta, 'phi': phi, 'Mx': Mx, 'My': My},
                    'outputs': {'A': A_matlab, 'dAtheta': dAtheta_matlab, 'dAphi': dAphi_matlab}
                }
                
            except Exception as e3:
                print(f"Error extracting data (method 3): {e3}")
                raise Exception(f"Could not extract data from MATLAB structure. Tried multiple methods.")

def run_test_case(case_data):
    """Run a single test case"""
    # Extract data using robust method
    extracted_data = extract_matlab_data(case_data)
    
    # Get inputs
    theta = extracted_data['inputs']['theta']
    phi = extracted_data['inputs']['phi']
    Mx = extracted_data['inputs']['Mx']
    My = extracted_data['inputs']['My']
    
    print(f"  Inputs: theta={theta}, phi={phi}, Mx={Mx}, My={My}")
    
    # Run Python implementation
    A_py, dAtheta_py, dAphi_py = construct_steer_matrix_and_derivative_steer_matrix(
        theta, phi, Mx, My
    )
    
    # Get MATLAB results
    A_matlab = extracted_data['outputs']['A']
    dAtheta_matlab = extracted_data['outputs']['dAtheta']
    dAphi_matlab = extracted_data['outputs']['dAphi']
    
    print(f"  Python A shape: {A_py.shape}, MATLAB A shape: {A_matlab.shape}")
    print(f"  Python dAtheta shape: {dAtheta_py.shape}, MATLAB dAtheta shape: {dAtheta_matlab.shape}")
    print(f"  Python dAphi shape: {dAphi_py.shape}, MATLAB dAphi shape: {dAphi_matlab.shape}")
    
    # Compare results
    A_match, A_diff, A_mean = compare_complex_matrices(A_py, A_matlab)
    dAtheta_match, dAtheta_diff, dAtheta_mean = compare_complex_matrices(dAtheta_py, dAtheta_matlab)
    dAphi_match, dAphi_diff, dAphi_mean = compare_complex_matrices(dAphi_py, dAphi_matlab)
    
    return {
        'python_results': (A_py, dAtheta_py, dAphi_py),
        'matlab_results': (A_matlab, dAtheta_matlab, dAphi_matlab),
        'A_match': A_match,
        'A_diff': A_diff,
        'A_mean': A_mean,
        'dAtheta_match': dAtheta_match,
        'dAtheta_diff': dAtheta_diff,
        'dAtheta_mean': dAtheta_mean,
        'dAphi_match': dAphi_match,
        'dAphi_diff': dAphi_diff,
        'dAphi_mean': dAphi_mean,
        'inputs': {
            'theta': theta,
            'phi': phi,
            'Mx': Mx,
            'My': My
        }
    }

def visualize_comparison(case_num, python_results, matlab_results, save_plots=True):
    """Create visualization comparing Python and MATLAB results"""
    python_A, python_dAtheta, python_dAphi = python_results
    matlab_A, matlab_dAtheta, matlab_dAphi = matlab_results
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle(f'Test Case {case_num} Comparison: Python vs MATLAB vs Difference')
    
    # Plot A
    im1 = axes[0,0].imshow(np.abs(python_A), cmap='viridis')
    axes[0,0].set_title('Python: |A|')
    plt.colorbar(im1, ax=axes[0,0])
    
    im2 = axes[0,1].imshow(np.abs(matlab_A), cmap='viridis')
    axes[0,1].set_title('MATLAB: |A|')
    plt.colorbar(im2, ax=axes[0,1])
    
    diff_A = np.abs(python_A - matlab_A)
    im3 = axes[0,2].imshow(diff_A, cmap='hot')
    axes[0,2].set_title(f'|Difference A| (max: {np.max(diff_A):.2e})')
    plt.colorbar(im3, ax=axes[0,2])
    
    # Plot dAtheta
    im4 = axes[1,0].imshow(np.abs(python_dAtheta), cmap='viridis')
    axes[1,0].set_title('Python: |dAtheta|')
    plt.colorbar(im4, ax=axes[1,0])
    
    im5 = axes[1,1].imshow(np.abs(matlab_dAtheta), cmap='viridis')
    axes[1,1].set_title('MATLAB: |dAtheta|')
    plt.colorbar(im5, ax=axes[1,1])
    
    diff_dAtheta = np.abs(python_dAtheta - matlab_dAtheta)
    im6 = axes[1,2].imshow(diff_dAtheta, cmap='hot')
    axes[1,2].set_title(f'|Difference dAtheta| (max: {np.max(diff_dAtheta):.2e})')
    plt.colorbar(im6, ax=axes[1,2])
    
    # Plot dAphi
    im7 = axes[2,0].imshow(np.abs(python_dAphi), cmap='viridis')
    axes[2,0].set_title('Python: |dAphi|')
    plt.colorbar(im7, ax=axes[2,0])
    
    im8 = axes[2,1].imshow(np.abs(matlab_dAphi), cmap='viridis')
    axes[2,1].set_title('MATLAB: |dAphi|')
    plt.colorbar(im8, ax=axes[2,1])
    
    diff_dAphi = np.abs(python_dAphi - matlab_dAphi)
    im9 = axes[2,2].imshow(diff_dAphi, cmap='hot')
    axes[2,2].set_title(f'|Difference dAphi| (max: {np.max(diff_dAphi):.2e})')
    plt.colorbar(im9, ax=axes[2,2])
    
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f'steer_matrix_comparison_case{case_num}.png', dpi=150, bbox_inches='tight')
        print(f"  Visualization saved as steer_matrix_comparison_case{case_num}.png")
    
    plt.close()

def main():
    """Main test function"""
    print("Testing construct_steer_matrix_and_derivative_steer_matrix")
    print("=" * 60)
    
    # Check if MATLAB test data exists
    matlab_data_file = '../../../OBS_for_CRLB_ISAC/tests/test_steer_matrix_data.mat'
    if not os.path.exists(matlab_data_file):
        print(f"Error: MATLAB test data not found at {matlab_data_file}")
        print("Please run the MATLAB script test_steer_matrix_generate_data.m first")
        return False
    
    # Load MATLAB test data
    try:
        matlab_data = loadmat(matlab_data_file)
        test_data = matlab_data['test_data']
        print("MATLAB test data loaded successfully")
        print(f"Test data type: {type(test_data)}")
        print(f"Test data shape: {test_data.shape}")
        if hasattr(test_data, 'dtype'):
            print(f"Test data dtype: {test_data.dtype}")
            if hasattr(test_data.dtype, 'names'):
                print(f"Test data field names: {test_data.dtype.names}")
    except Exception as e:
        print(f"Error loading MATLAB test data: {e}")
        return False
    
    # Run tests for each case
    all_passed = True
    results_summary = []
    
    # Try to access test cases
    try:
        # Method 1: Access as structured array
        if hasattr(test_data.dtype, 'names') and test_data.dtype.names:
            case_names = [name for name in test_data.dtype.names if name.startswith('case')]
        else:
            # Method 2: Assume cases 1, 2, 3 exist
            case_names = ['case1', 'case2', 'case3']
        
        for case_name in case_names:
            case_num = int(case_name.replace('case', ''))
            
            try:
                # Try different ways to access the case data
                if hasattr(test_data.dtype, 'names') and case_name in test_data.dtype.names:
                    case_data = test_data[case_name][0, 0]
                else:
                    print(f"Warning: {case_name} not found in test data")
                    continue
                
                print(f"\nTest Case {case_num}:")
                print("-" * 20)
                
                # Run test case
                results = run_test_case(case_data)
                
                # Print results
                print(f"  Array size: {results['inputs']['Mx']}x{results['inputs']['My']}")
                print(f"  Number of angles: {len(results['inputs']['theta'])}")
                print(f"  A matrices match: {results['A_match']} (max diff: {results['A_diff']:.2e}, mean diff: {results['A_mean']:.2e})")
                print(f"  dAtheta matrices match: {results['dAtheta_match']} (max diff: {results['dAtheta_diff']:.2e}, mean diff: {results['dAtheta_mean']:.2e})")
                print(f"  dAphi matrices match: {results['dAphi_match']} (max diff: {results['dAphi_diff']:.2e}, mean diff: {results['dAphi_mean']:.2e})")
                
                # Check if all matrices match
                case_passed = results['A_match'] and results['dAtheta_match'] and results['dAphi_match']
                print(f"  Case {case_num} overall: {'PASS' if case_passed else 'FAIL'}")
                
                if not case_passed:
                    all_passed = False
                
                # Create visualization
                visualize_comparison(case_num, results['python_results'], results['matlab_results'])
                
                # Store summary
                results_summary.append({
                    'case': case_num,
                    'passed': case_passed,
                    'max_errors': {
                        'A': results['A_diff'],
                        'dAtheta': results['dAtheta_diff'],
                        'dAphi': results['dAphi_diff']
                    },
                    'mean_errors': {
                        'A': results['A_mean'],
                        'dAtheta': results['dAtheta_mean'],
                        'dAphi': results['dAphi_mean']
                    }
                })
                
            except Exception as e:
                print(f"Error processing {case_name}: {e}")
                import traceback
                traceback.print_exc()
                all_passed = False
    
    except Exception as e:
        print(f"Error accessing test cases: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Overall summary
    print(f"\n{'='*60}")
    print("OVERALL TEST SUMMARY")
    print(f"{'='*60}")
    print(f"All tests passed: {'YES' if all_passed else 'NO'}")
    
    if all_passed:
        print("✅ The Python implementation matches the MATLAB implementation!")
    else:
        print("❌ There are differences between Python and MATLAB implementations.")
        print("Check the visualizations and error details above.")
    
    # Save detailed results
    try:
        savemat('steer_matrix_test_results.mat', {
            'results_summary': results_summary,
            'all_passed': all_passed
        })
        print(f"\nDetailed results saved to steer_matrix_test_results.mat")
    except Exception as e:
        print(f"Warning: Could not save results: {e}")
    
    print("Visualizations saved as steer_matrix_comparison_case*.png")
    
    return all_passed

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
