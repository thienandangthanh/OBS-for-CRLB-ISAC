#!/usr/bin/env python3
"""
Modify .mat file metadata for testing purposes
Takes input and output file paths as command line arguments
"""

import scipy.io
from scipy.io import savemat
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import argparse

def modify_mat_metadata(input_file, output_file):
    """
    Load a .mat file and save it with potentially different metadata

    Args:
        input_file (str): Path to input .mat file
        output_file (str): Path to output .mat file
    """
    # Check if input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Load the MATLAB data file
    print(f"Loading data from: {input_file}")
    data = scipy.io.loadmat(input_file)

    # Print the keys to see what variables are available
    print("Available variables in the .mat file:")
    print(list(data.keys()))

    # Remove metadata and keep only data variables
    clean_data = {k: v for k, v in data.items() if not k.startswith('__')}

    print(f"Data variables to save: {list(clean_data.keys())}")

    # Print info about each variable
    for var_name, var_data in clean_data.items():
        if hasattr(var_data, 'shape'):
            print(f"{var_name} shape: {var_data.shape}, type: {type(var_data)}")
        else:
            print(f"{var_name} type: {type(var_data)}")

    # Save to output file - this will create new metadata
    print(f"Saving modified data to: {output_file}")
    savemat(output_file, clean_data)

    print(f"✓ Modified data saved successfully to: {output_file}")
    print("  Note: New metadata will be generated automatically by scipy.io.savemat")

def main():
    parser = argparse.ArgumentParser(
        description='Modify .mat file metadata for testing purposes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python mod-mat-metadata.py input.mat output.mat
        """
    )

    parser.add_argument('input', help='Path to input .mat file')
    parser.add_argument('output', help='Path to output .mat file')

    args = parser.parse_args()

    try:
        modify_mat_metadata(args.input, args.output)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
