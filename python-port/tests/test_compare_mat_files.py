#!/usr/bin/env python3
"""
Test script for compare_mat_files.py
Demonstrates that metadata differences don't affect data comparison results
"""

import os
import sys
import subprocess

# Add parent directory to path to import compare_mat_files
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'helpers'))
from compare_mat_files import compare_mat_files

def test_metadata_differences():
    """
    Test that compare_mat_files.py correctly handles metadata differences
    """
    print("Running test: Metadata differences should not affect data comparison")
    print("=" * 70)

    # Step 1: Generate modified metadata file
    print("Step 1: Generating modified metadata file...")
    mod_script_path = os.path.join(os.path.dirname(__file__), 'mod-mat-metadata.py')

    # Define input and output file paths
    original_file = os.path.join(os.path.dirname(__file__), '..', 'Fig1_convergence', 'data_convergence.mat')
    modified_file = os.path.join(os.path.dirname(__file__), 'mod_data_convergence.mat')

    # Check if original file exists before proceeding
    if not os.path.exists(original_file):
        print(f"Error: Original file not found: {original_file}")
        return False

    try:
        # Run the mod-mat-metadata.py script with input and output arguments
        result = subprocess.run([sys.executable, mod_script_path, original_file, modified_file], 
                              capture_output=True, text=True, cwd=os.path.dirname(mod_script_path))
        if result.returncode != 0:
            print(f"Error running mod-mat-metadata.py: {result.stderr}")
            return False
        print("✓ Modified metadata file generated successfully")
        print(result.stdout)
    except Exception as e:
        print(f"Error running mod-mat-metadata.py: {e}")
        return False

    # Step 2: Compare original and modified files
    print("\nStep 2: Comparing original and modified files...")

    if not os.path.exists(modified_file):
        print(f"Error: Modified file not found: {modified_file}")
        return False

    try:
        # Run comparison with verbose output
        summary = compare_mat_files(
            original_file,
            modified_file,
            tolerance_abs=1e-12,
            tolerance_rel=1e-9,
            verbose=True
        )

        # Check results
        if summary['success']:
            print("\n" + "=" * 70)
            print("TEST RESULTS:")
            print("=" * 70)

            # Check data comparison results
            data_failed = summary['failed']
            metadata_different = len(summary['metadata_comparison']['differences'])

            print(f"Data comparisons failed: {data_failed}")
            print(f"Metadata fields different: {metadata_different}")

            if data_failed == 0 and metadata_different > 0:
                print("✅ TEST PASSED: Data is identical despite metadata differences")
                print("   This confirms that compare_mat_files.py correctly separates")
                print("   metadata differences from data differences.")
                return True
            elif data_failed == 0 and metadata_different == 0:
                print("⚠️  WARNING: No metadata differences detected")
                print("   The modified file may be identical to the original")
                return True
            else:
                print("❌ TEST FAILED: Data comparison failed unexpectedly")
                return False
        else:
            print("❌ TEST FAILED: File comparison failed")
            return False

    except Exception as e:
        print(f"Error during comparison: {e}")
        return False

def cleanup():
    """Clean up generated test files"""
    test_file = os.path.join(os.path.dirname(__file__), 'mod_data_convergence.mat')
    if os.path.exists(test_file):
        os.remove(test_file)
        print(f"Cleaned up: {test_file}")

if __name__ == "__main__":
    try:
        success = test_metadata_differences()

        if success:
            print("\n🎉 All tests passed!")
        else:
            print("\n❌ Tests failed!")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\nTest interrupted by user")
    finally:
        # Always clean up
        print("\nCleaning up...")
        cleanup()
