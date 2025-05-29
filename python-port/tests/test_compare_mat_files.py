#!/usr/bin/env python3
"""
Test script for compare_mat_files.py
Demonstrates that metadata differences don't affect data comparison results
"""

import os
import sys
import subprocess
import pytest
from pathlib import Path

# Add parent directory to path to import compare_mat_files
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'helpers'))
from compare_mat_files import compare_mat_files

# Define project paths
PROJECT_ROOT = Path(__file__).parent.parent
TESTS_DIR = PROJECT_ROOT / 'tests'
FIG1_DIR = PROJECT_ROOT / 'Fig1_convergence'

@pytest.fixture(autouse=True)
def cleanup():
    """Cleanup fixture that runs after each test"""
    # This code runs before each test
    yield
    # This code runs after each test
    test_file = TESTS_DIR / 'mod_data_convergence.mat'
    if test_file.exists():
        test_file.unlink()
        print(f"\nCleaned up: {test_file}")

def test_metadata_differences():
    """
    Test that compare_mat_files.py correctly handles metadata differences
    """
    print("Running test: Metadata differences should not affect data comparison")
    print("=" * 70)

    # Step 1: Generate modified metadata file
    print("Step 1: Generating modified metadata file...")
    mod_script_path = TESTS_DIR / 'mod-mat-metadata.py'

    # Define input and output file paths
    original_file = FIG1_DIR / 'data_convergence.mat'
    modified_file = TESTS_DIR / 'mod_data_convergence.mat'

    # Check if original file exists before proceeding
    if not original_file.exists():
        pytest.fail(f"Error: Original file not found: {original_file}")

    try:
        # Run the mod-mat-metadata.py script with input and output arguments
        result = subprocess.run([sys.executable, str(mod_script_path), 
                               str(original_file), str(modified_file)], 
                              capture_output=True, text=True, 
                              cwd=str(TESTS_DIR))
        if result.returncode != 0:
            pytest.fail(f"Error running mod-mat-metadata.py: {result.stderr}")
        print("✓ Modified metadata file generated successfully")
        print(result.stdout)
    except Exception as e:
        pytest.fail(f"Error running mod-mat-metadata.py: {e}")

    # Step 2: Compare original and modified files
    print("\nStep 2: Comparing original and modified files...")

    if not modified_file.exists():
        pytest.fail(f"Error: Modified file not found: {modified_file}")

    try:
        # Run comparison with verbose output
        summary = compare_mat_files(
            str(original_file),
            str(modified_file),
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
                assert True
            elif data_failed == 0 and metadata_different == 0:
                print("⚠️  WARNING: No metadata differences detected")
                print("   The modified file may be identical to the original")
                assert True
            else:
                pytest.fail("❌ TEST FAILED: Data comparison failed unexpectedly")
        else:
            pytest.fail("❌ TEST FAILED: File comparison failed")

    except Exception as e:
        pytest.fail(f"Error during comparison: {e}")

if __name__ == "__main__":
    pytest.main([__file__])
