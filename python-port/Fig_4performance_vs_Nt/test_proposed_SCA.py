#!/usr/bin/env python3
"""
Test script for Figure 4 proposed_SCA.py implementation.

This script runs a quick test with reduced parameters to validate 
the implementation before running full simulations.
"""

import sys
import os

# Add the current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_proposed_SCA():
    """Run a quick test of the proposed SCA implementation."""
    
    print("Testing Figure 4 proposed_SCA.py implementation...")
    print("=" * 60)
    
    # Import the main function from proposed_SCA
    from proposed_SCA import main
    
    # Override sys.argv to provide test arguments
    original_argv = sys.argv
    try:
        sys.argv = [
            'proposed_SCA.py',
            '--i-out', '3',          # Reduced channel realizations 
            '--save-data',           # Save test data
            '--save-plots',          # Save test plots
            '--output-dir', './test_results',
            '--tolerance', '1e-4',   # Relaxed tolerance for faster convergence
            '--max-iterations', '100' # Reduced iterations for faster test
        ]
        
        # Create output directory
        os.makedirs('./test_results', exist_ok=True)
        
        # Run the main function
        print("Running proposed_SCA with test parameters...")
        SR_all, CRB_all, Time_all, Obj_all, Con = main()
        
        # Validate results
        print("\nValidation Results:")
        print("-" * 30)
        
        # Check dimensions
        expected_shape = (3, 5)  # (I_out=3, I_in=5)
        assert SR_all.shape == expected_shape, f"SR_all shape mismatch: {SR_all.shape} vs {expected_shape}"
        assert CRB_all.shape == expected_shape, f"CRB_all shape mismatch: {CRB_all.shape} vs {expected_shape}"
        assert Time_all.shape == expected_shape, f"Time_all shape mismatch: {Time_all.shape} vs {expected_shape}"
        assert Obj_all.shape == expected_shape, f"Obj_all shape mismatch: {Obj_all.shape} vs {expected_shape}"
        print("✓ Array dimensions are correct")
        
        # Check for finite values
        assert all([
            SR_all.max() < float('inf'),
            CRB_all.max() < float('inf'),
            Time_all.max() < float('inf'),
            Obj_all.max() < float('inf')
        ]), "Some results contain infinite values"
        print("✓ All values are finite")
        
        # Check for reasonable ranges
        assert SR_all.min() >= 0, f"Sum rate cannot be negative: {SR_all.min()}"
        assert CRB_all.min() > 0, f"CRB trace must be positive: {CRB_all.min()}"
        assert Time_all.min() > 0, f"Computation time must be positive: {Time_all.min()}"
        print("✓ Value ranges are reasonable")
        
        # Check convergence data
        assert len(Con) > 0, "No convergence data recorded"
        assert len(Con[0]) == 3, f"Convergence data format incorrect: {len(Con[0])} columns"
        print("✓ Convergence data is valid")
        
        # Print summary statistics
        print("\nSummary Statistics:")
        print(f"  Sum Rate: {SR_all.mean():.4f} ± {SR_all.std():.4f}")
        print(f"  CRB Trace: {CRB_all.mean():.4e} ± {CRB_all.std():.4e}")
        print(f"  Objective: {Obj_all.mean():.4f} ± {Obj_all.std():.4f}")
        print(f"  Time per scenario: {Time_all.mean():.3f}s ± {Time_all.std():.3f}s")
        print(f"  Convergence iterations: {len(Con)} (last scenario)")
        
        print("\n✅ Test completed successfully!")
        print("You can now run the full simulation with:")
        print("python proposed_SCA.py --save-data --save-plots")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        # Restore original sys.argv
        sys.argv = original_argv


if __name__ == "__main__":
    success = test_proposed_SCA()
    exit(0 if success else 1) 