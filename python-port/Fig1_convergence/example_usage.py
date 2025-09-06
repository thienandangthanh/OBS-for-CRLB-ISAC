#!/usr/bin/env python3
"""
Example usage of the refactored propose_SCA.py script

This script demonstrates different ways to run the convergence analysis
with various configurations.
"""

import subprocess
import os

def run_example(description, args):
    """Run an example with the given arguments"""
    print(f"\n{'='*60}")
    print(f"Example: {description}")
    print(f"{'='*60}")
    print(f"Command: python propose_SCA.py {' '.join(args)}")
    print()

    # Run the command
    result = subprocess.run(['python', 'propose_SCA.py'] + args,
                            capture_output=True, text=True)

    if result.returncode == 0:
        print("✓ Success!")
        # Show only the last few lines of output
        lines = result.stdout.strip().split('\n')
        if len(lines) > 10:
            print("Output (last 10 lines):")
            for line in lines[-10:]:
                print(f"  {line}")
        else:
            print("Output:")
            print(result.stdout)
    else:
        print("✗ Failed!")
        print("Error:")
        print(result.stderr)

def main():
    print("OBS-for-CRLB-ISAC Figure 1 Convergence Analysis Examples")
    print("======================================================")

    # Example 1: Default configuration with saved plots
    run_example(
        "Default configuration with saved plots",
        ['--save-plots', '--output-dir', './example_output1', '--max-iterations', '10']
    )

    # Example 2: Modified antenna configuration
    run_example(
        "Modified antenna configuration (6x6 transmit)",
        ['--nth', '6', '--ntv', '6', '--save-plots', '--output-dir', './example_output2', '--max-iterations', '10']
    )

    # Example 3: Different power settings
    run_example(
        "Higher transmit power and different noise levels",
        ['--pt-dbm', '15', '--noise-c-dbm', '-5', '--noise-s-dbm', '-5', 
         '--save-plots', '--output-dir', './example_output3', '--max-iterations', '10']
    )

    # Example 4: More users and targets
    run_example(
        "More communication users and sensing targets",
        ['--k', '6', '--m', '3', '--save-plots', '--output-dir', './example_output4', '--max-iterations', '10']
    )

    # Example 5: Different convergence settings
    run_example(
        "Tighter convergence tolerance",
        ['--tolerance', '1e-6', '--max-iterations', '20', '--save-plots', '--output-dir', './example_output5']
    )

    # Example 6: Complete workflow with combined plot
    print(f"\n{'='*60}")
    print("Example: Complete workflow with combined convergence plot")
    print(f"{'='*60}")

    # First run the main script
    print("Step 1: Running main convergence analysis...")
    result1 = subprocess.run([
        'python', 'propose_SCA.py', 
        '--save-plots', '--output-dir', './example_output6', '--max-iterations', '10'
    ], capture_output=True, text=True)

    if result1.returncode == 0:
        print("✓ Main analysis completed!")

        # Then generate the combined plot
        print("Step 2: Generating combined convergence plot...")
        result2 = subprocess.run([
            'python', 'plotfigure1_convergence.py',
            '--input-file', './example_output6/data_convergence.mat',
            '--save-only', '--output-dir', './example_output6'
        ], capture_output=True, text=True)

        if result2.returncode == 0:
            print("✓ Combined plot generated!")
            print("Output files in ./example_output6/:")
            print("  - data_convergence.mat")
            print("  - fig1_sum_rate_vs_iteration.png")
            print("  - fig1_inverse_fim_trace_vs_iteration.png") 
            print("  - fig1_objective_value_vs_iteration.png")
            print("  - fig1_convergence_combined.png")
        else:
            print("✗ Combined plot generation failed!")
            print("Error:", result2.stderr)
    else:
        print("✗ Main analysis failed!")
        print("Error:", result1.stderr)

    print(f"\n{'='*60}")
    print("All examples completed!")
    print("Check the example_output* directories for results and plots.")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
