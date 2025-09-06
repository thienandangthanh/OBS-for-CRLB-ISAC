import numpy as np
import pandas as pd
from scipy.io import loadmat
import os

def complex_to_string(x):
    """Convert complex number to readable string format"""
    if np.imag(x) == 0:
        return f"{np.real(x):.6f}"
    elif np.real(x) == 0:
        return f"{np.imag(x):.6f}j"
    else:
        return f"{np.real(x):.6f}{'+' if np.imag(x) >= 0 else ''}{np.imag(x):.6f}j"

def handle_multidimensional_array(var_data, var_name, writer):
    """Handle arrays with more than 2 dimensions by flattening or slicing"""
    shape = var_data.shape
    print(f"  Handling {len(shape)}D array with shape {shape}")

    if len(shape) == 3:
        # For 3D arrays, create separate sheets for each "page" along the 3rd dimension
        for i in range(shape[2]):
            slice_data = var_data[:, :, i]

            # Convert complex numbers to string if needed
            if np.iscomplexobj(slice_data):
                slice_data = np.vectorize(complex_to_string)(slice_data)

            # Create DataFrame
            df = pd.DataFrame(slice_data)
            sheet_name = f"{var_name}_dim{i}"
            print(f"    Creating sheet: {sheet_name}")
            df.to_excel(writer, sheet_name=sheet_name)

    elif len(shape) == 4:
        # For 4D arrays, create sheets for each combination of the last two dimensions
        for i in range(shape[2]):
            for j in range(shape[3]):
                slice_data = var_data[:, :, i, j]

                # Convert complex numbers to string if needed
                if np.iscomplexobj(slice_data):
                    slice_data = np.vectorize(complex_to_string)(slice_data)

                # Create DataFrame
                df = pd.DataFrame(slice_data)
                sheet_name = f"{var_name}_dim{i}_{j}"
                print(f"    Creating sheet: {sheet_name}")
                df.to_excel(writer, sheet_name=sheet_name)

    else:
        # For higher dimensions, flatten to 2D by reshaping
        print(f"    Flattening {len(shape)}D array to 2D")

        # Reshape to 2D by combining all dimensions except the first
        reshaped_data = var_data.reshape(shape[0], -1)

        # Convert complex numbers to string if needed
        if np.iscomplexobj(reshaped_data):
            reshaped_data = np.vectorize(complex_to_string)(reshaped_data)

        # Create DataFrame
        df = pd.DataFrame(reshaped_data)
        sheet_name = f"{var_name}_flattened"
        print(f"    Creating sheet: {sheet_name}")
        df.to_excel(writer, sheet_name=sheet_name)

def mat_to_excel(mat_file, excel_file=None):
    """
    Convert .mat file to .xlsx file

    Args:
        mat_file (str): Path to .mat file
        excel_file (str, optional): Path to output .xlsx file. If None, uses same name as mat_file
    """
    # Load the .mat file
    print(f"Loading {mat_file}...")
    data = loadmat(mat_file)

    # If excel_file not specified, use same name as mat_file but with .xlsx extension
    if excel_file is None:
        excel_file = os.path.splitext(mat_file)[0] + '.xlsx'

    print(f"Converting to {excel_file}...")

    # Track if any sheets were created
    sheets_created = 0

    # Create Excel writer
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        # Process each variable in the .mat file
        for var_name, var_data in data.items():
            # Skip MATLAB's internal variables
            if var_name.startswith('__'):
                continue

            print(f"Processing variable: {var_name}")

            # Handle different data types
            if isinstance(var_data, np.ndarray):
                # Check dimensions
                if var_data.ndim == 0:
                    # Scalar values
                    if np.iscomplexobj(var_data):
                        var_data = complex_to_string(var_data.item())
                    df = pd.DataFrame([[var_data.item()]], columns=[var_name])
                    print(f"  Creating sheet: {var_name}")
                    df.to_excel(writer, sheet_name=var_name, index=False)
                    sheets_created += 1

                elif var_data.ndim == 1:
                    # 1D arrays - convert to column vector
                    if np.iscomplexobj(var_data):
                        var_data = np.vectorize(complex_to_string)(var_data)
                    df = pd.DataFrame(var_data.reshape(-1, 1), columns=[var_name])
                    print(f"  Creating sheet: {var_name}")
                    df.to_excel(writer, sheet_name=var_name)
                    sheets_created += 1

                elif var_data.ndim == 2:
                    # For cell arrays in MATLAB (stored as object arrays in Python)
                    if var_data.dtype == np.dtype('O'):
                        for i, cell_data in enumerate(var_data.flatten()):
                            if hasattr(cell_data, 'size') and cell_data.size > 0:  # Only process non-empty cells
                                # Convert complex numbers to string if needed
                                if np.iscomplexobj(cell_data):
                                    cell_data = np.vectorize(complex_to_string)(cell_data)

                                # Create DataFrame
                                df = pd.DataFrame(cell_data)
                                sheet_name = f"{var_name}_{i}"
                                print(f"  Creating sheet: {sheet_name}")
                                df.to_excel(writer, sheet_name=sheet_name)
                                sheets_created += 1

                        # For regular 2D arrays
                    else:
                        # Convert complex numbers to string if needed
                        if np.iscomplexobj(var_data):
                            var_data = np.vectorize(complex_to_string)(var_data)

                        # Create DataFrame
                        df = pd.DataFrame(var_data)
                        print(f"  Creating sheet: {var_name}")
                        df.to_excel(writer, sheet_name=var_name)
                        sheets_created += 1

                else:
                    # Handle multidimensional arrays (3D, 4D, etc.)
                    handle_multidimensional_array(var_data, var_name, writer)
                    sheets_created += 1

            else:
                # Handle other data types (strings, etc.)
                df = pd.DataFrame([[str(var_data)]], columns=[var_name])
                print(f"  Creating sheet: {var_name}")
                df.to_excel(writer, sheet_name=var_name, index=False)
                sheets_created += 1

    if sheets_created == 0:
        print("Warning: No valid data found to convert!")
        # Create a dummy sheet to avoid "At least one sheet must be visible" error
        df = pd.DataFrame([["No valid data found"]], columns=["Info"])
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name="NoData", index=False)
        return

    print("Conversion completed successfully!")

    # Print basic information about the converted data
    print("\nData Summary:")
    try:
        excel_data = pd.read_excel(excel_file, sheet_name=None)
        for sheet_name, df in list(excel_data.items())[:5]:  # Show first 5 sheets only
            print(f"\nSheet: {sheet_name}")
            print(f"Shape: {df.shape}")
            print("First few rows:")
            print(df.head())

        if len(excel_data) > 5:
            print(f"\n... and {len(excel_data) - 5} more sheets")

    except Exception as e:
        print(f"Could not read back the Excel file for summary: {e}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Convert MATLAB .mat file to Excel .xlsx file')
    parser.add_argument('mat_file', help='Path to input .mat file')
    parser.add_argument('--output', '-o', help='Path to output .xlsx file (optional)')

    args = parser.parse_args()

    mat_to_excel(args.mat_file, args.output)
