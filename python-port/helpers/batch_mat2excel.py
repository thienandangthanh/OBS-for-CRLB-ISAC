#!/usr/bin/env python3
"""
Batch converter for .mat files to .xlsx files
Finds all .mat files in a directory (with exclusions) and converts them using mat2excel.py
"""

import os
import argparse
import sys
from pathlib import Path

# Import the mat2excel function from the same directory
from mat2excel import mat_to_excel

def find_mat_files(include_dir, exclude_dirs=None):
    """
    Find all .mat files in include_dir, excluding files in exclude_dirs

    Args:
        include_dir (str): Directory to search for .mat files
        exclude_dirs (list): List of directories to exclude from search

    Returns:
        list: List of paths to .mat files
    """
    if exclude_dirs is None:
        exclude_dirs = []

    # Convert to Path objects and resolve to absolute paths
    include_path = Path(include_dir).resolve()
    exclude_paths = [Path(exclude_dir).resolve() for exclude_dir in exclude_dirs]

    print(f"Searching for .mat files in: {include_path}")
    if exclude_paths:
        print(f"Excluding directories: {[str(p) for p in exclude_paths]}")

    mat_files = []

    # Walk through all subdirectories
    for root, dirs, files in os.walk(include_path):
        root_path = Path(root).resolve()

        # Check if current directory should be excluded
        skip_directory = False
        for exclude_path in exclude_paths:
            try:
                # Check if current directory is within any excluded directory
                root_path.relative_to(exclude_path)
                skip_directory = True
                print(f"Skipping excluded directory: {root_path}")
                break
            except ValueError:
                # Not within this excluded directory, continue checking
                continue

        if skip_directory:
            # Clear dirs to prevent os.walk from descending into subdirectories
            dirs.clear()
            continue

        # Find .mat files in current directory
        for file in files:
            if file.lower().endswith('.mat'):
                mat_file_path = os.path.join(root, file)
                mat_files.append(mat_file_path)
                print(f"Found .mat file: {mat_file_path}")

    return mat_files

def batch_convert_mat_files(include_dir, exclude_dirs=None):
    """
    Convert all .mat files to .xlsx files

    Args:
        include_dir (str): Directory to search for .mat files
        exclude_dirs (list): List of directories to exclude from search
    """
    # Find all .mat files
    mat_files = find_mat_files(include_dir, exclude_dirs)

    if not mat_files:
        print("No .mat files found!")
        return

    print(f"\nFound {len(mat_files)} .mat file(s) to convert:")
    for mat_file in mat_files:
        print(f"  - {mat_file}")

    # Convert each .mat file
    successful_conversions = 0
    failed_conversions = 0

    for i, mat_file in enumerate(mat_files, 1):
        print(f"\n{'='*60}")
        print(f"Converting file {i}/{len(mat_files)}: {mat_file}")
        print(f"{'='*60}")

        try:
            # The mat_to_excel function will create the .xlsx file in the same directory
            # as the .mat file with the same name but .xlsx extension
            mat_to_excel(mat_file)
            successful_conversions += 1
            print(f"✓ Successfully converted: {mat_file}")

        except Exception as e:
            failed_conversions += 1
            print(f"✗ Failed to convert {mat_file}: {str(e)}")

    # Summary
    print(f"\n{'='*60}")
    print("CONVERSION SUMMARY")
    print(f"{'='*60}")
    print(f"Total files found: {len(mat_files)}")
    print(f"Successful conversions: {successful_conversions}")
    print(f"Failed conversions: {failed_conversions}")

    if failed_conversions > 0:
        print(f"\n⚠️  {failed_conversions} file(s) failed to convert. Check the error messages above.")
    else:
        print(f"\n🎉 All files converted successfully!")

def main():
    parser = argparse.ArgumentParser(
        description='Batch convert .mat files to .xlsx files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert all .mat files in current directory
  python batch_mat2excel.py --include .

  # Convert all .mat files in data/ directory, excluding temp/ subdirectory
  python batch_mat2excel.py --include data --exclude temp

  # Convert with multiple exclusions
  python batch_mat2excel.py --include /path/to/data --exclude temp backup old_files
        """
    )

    parser.add_argument(
        '--include', 
        required=True,
        help='Directory to search for .mat files'
    )

    parser.add_argument(
        '--exclude',
        nargs='*',
        default=[],
        help='Directory(ies) to exclude from search (can specify multiple)'
    )

    args = parser.parse_args()

    # Validate include directory exists
    if not os.path.exists(args.include):
        print(f"Error: Include directory '{args.include}' does not exist!")
        sys.exit(1)

    if not os.path.isdir(args.include):
        print(f"Error: Include path '{args.include}' is not a directory!")
        sys.exit(1)

    # Validate exclude directories exist (if specified)
    for exclude_dir in args.exclude:
        if not os.path.exists(exclude_dir):
            print(f"Warning: Exclude directory '{exclude_dir}' does not exist, ignoring...")
            args.exclude.remove(exclude_dir)
        elif not os.path.isdir(exclude_dir):
            print(f"Warning: Exclude path '{exclude_dir}' is not a directory, ignoring...")
            args.exclude.remove(exclude_dir)

    # Perform batch conversion
    batch_convert_mat_files(args.include, args.exclude if args.exclude else None)

if __name__ == "__main__":
    main()
