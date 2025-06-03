#!/usr/bin/env python3
"""
Example usage of the refactored monte_carlo_verification_python.py script

This script demonstrates different ways to run Monte Carlo verification
with various configurations.
"""

import subprocess
import os

def run_mc_example(description, args):
    """Run a Monte Carlo example with the given arguments"""
    print(f"\n{'='*80}")
    print(f"Example: {description}")
    print(f"{'='*80}")
    print(f"Command: python monte_carlo_verification_python.py {' '.join(args)}")
    print()
    
    # Run the command
    result = subprocess.run(['python', 'monte_carlo_verification_python.py'] + args, 
                          capture_output=True, text=True)
    
    if result.returncode == 0:
        print("✓ Success!")
        # Show only the last few lines of output
        lines = result.stdout.strip().split('\n')
        if len(lines) > 15:
            print("Output (last 15 lines):")
            for line in lines[-15:]:
                print(f"  {line}")
        else:
            print("Output:")
            print(result.stdout)
    else:
        print("✗ Failed!")
        print("Error:")
        print(result.stderr)

def main():
    print("Monte Carlo Verification Examples")
    print("==================================")
    print("These examples demonstrate the refactored Monte Carlo verification")
    print("script with different configurations using SimulationConfig.")
    
    # Example 1: Quick test with default configuration
    run_mc_example(
        "Quick test with 5 runs and limited iterations",
        ['--num-runs', '5', '--max-iterations', '20', '--output-dir', './mc_example1']
    )
    
    # Example 2: Different antenna configuration
    run_mc_example(
        "Modified antenna configuration (6x6 transmit)",
        ['--num-runs', '3', '--nth', '6', '--ntv', '6', '--max-iterations', '15', 
         '--output-dir', './mc_example2']
    )
    
    # Example 3: Different system parameters
    run_mc_example(
        "More users and targets with higher power",
        ['--num-runs', '4', '--k', '6', '--m', '3', '--pt-dbm', '15', 
         '--max-iterations', '25', '--output-dir', './mc_example3']
    )
    
    # Example 4: Different convergence settings
    run_mc_example(
        "Tighter convergence tolerance",
        ['--num-runs', '3', '--tolerance', '1e-6', '--max-iterations', '30', 
         '--output-dir', './mc_example4']
    )
    
    # Example 5: Save individual runs
    run_mc_example(
        "Save individual run data",
        ['--num-runs', '3', '--max-iterations', '15', '--save-individual', 
         '--output-dir', './mc_example5']
    )
    
    print(f"\n{'='*80}")
    print("All Monte Carlo examples completed!")
    print("Check the mc_example* directories for results and plots.")
    print(f"{'='*80}")
    
    # Summary of what was generated
    print("\nFiles generated in each example directory:")
    print("  - monte_carlo_averaged_results_python.mat")
    print("  - monte_carlo_averaged_convergence_python.png")
    print("  - monte_carlo_statistics_python.png")
    print("  - monte_carlo_individual_runs_python.mat (if --save-individual used)")

if __name__ == "__main__":
    main() 