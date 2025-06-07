# Figure 4: Performance vs Number of Transmit Antennas (Nt)

This directory contains the Python translation of the MATLAB Figure 4 analysis for the OBS-for-CRLB-ISAC system.

## Overview

Figure 4 analyzes how the system performance varies with the number of transmit antennas (Nt) while keeping other parameters fixed. Unlike Figures 2 and 3 which analyze trade-offs and sensing streams respectively, Figure 4 focuses on the scalability of the system with antenna array size.

## Key Differences Between SCA Implementations

The following table summarizes the key differences between the SCA algorithm implementations across Figures 2, 3, and 4:

| Parameter | Figure 2 SCA | Figure 3 SCA | Figure 4 SCA |
|-----------|--------------|--------------|--------------|
| **Primary Variable** | Trade-off weight (δc) | Number of sensing streams (Ns) | Number of antennas (Nt, Nr) |
| **M_max** | 2 | 3 | 2 |
| **m_target** | 2 | 3 | 2 |
| **K (users)** | 4 (fixed) | Loop: 1,2,3 | 4 (fixed) |
| **L (snapshots)** | 30 | 128 | 30 |
| **Antenna Config** | Fixed: Nth=4, Ntv=4 | Fixed: Nth=4, Ntv=4 | Variable: Nth=4, Ntv=2^(weight+1) |
| **Receive Config** | Fixed: Nrh=5, Nrv=4 | Fixed: Nrh=5, Nrv=4 | Variable: Nrh=5, Nrv=2^(weight+1) |
| **num_sensing_streams** | Nt (fixed) | weight-1 (0-8) | 3*M = 6 (fixed) |
| **delta_s** | 1 (fixed) | 1 (fixed) | 1 (fixed) |
| **delta_c** | 10^(delta_all(weight)) | 0.25 (fixed) | 0.25 (fixed) |
| **H_all dimensions** | (I_out, Nt, K) | (Nt, K_max, I_out) | (I_out, Nt, K) |
| **Wc initialization** | Random | delta_c*H/‖H‖ | delta_c*H/‖H‖ |
| **Main loop structure** | Single parameter sweep | User loop + parameter sweep | Single parameter sweep |
| **Tolerance** | 1e-4 | 1e-4 | 1e-5 |
| **Max iterations** | 2000 | 2000 | 4000 |
| **Variable range** | δc: -7 to 4.8 dB (60 points) | Ns: 0-8 streams | Nt: 8,16,32,64,128 antennas |

## Implementation Scope

This implementation includes:
- The main SCA algorithm with variable antenna arrays (proposed_SCA.py)
- Performance analysis versus number of transmit antennas
- Comparison with other algorithms (WMMSE-SDR, FP-SGDA, LD Algorithm)
- Data preservation from MATLAB implementation for algorithm comparison

## Key Features of Figure 4 Implementation

### Variable Antenna Configuration
- **Transmit**: Nth = 4, Ntv = 2^(weight+1), so Nt = 4 × 2^(weight+1)
- **Receive**: Nrh = 5, Nrv = 2^(weight+1), so Nr = 5 × 2^(weight+1)
- **Antenna values**: 8, 16, 32, 64, 128 (for weight = 1,2,3,4,5)

### Fixed System Parameters
- **Users**: K = 4 (fixed)
- **Targets**: M = 2 (fixed)
- **Sensing streams**: num_sensing_streams = 3*M = 6 (fixed)
- **Trade-off weights**: delta_c = 0.25, delta_s = 1 (fixed)

### Algorithm Configuration
- **Tolerance**: 1e-5 (tighter than Fig2/Fig3)
- **Max iterations**: 4000 (more than Fig2/Fig3)
- **Communication beamforming**: Wc = delta_c*H/‖H‖ (normalized per column)

## Usage

### Basic Usage
```bash
python proposed_SCA.py
```

### Quick Test (Reduced Parameters)
```bash
python proposed_SCA.py --i-out 10 --save-data --save-plots
```

### Custom Configuration
```bash
python proposed_SCA.py --tolerance 1e-6 --max-iterations 2000 --output-dir ./results
```

## Command Line Options

All standard simulation parameters are available via command line. Key options:

- `--i-out`: Number of channel realizations (default: 100)
- `--save-data`: Save results to data files
- `--save-plots`: Save plots as files instead of displaying
- `--output-dir`: Directory for output files (default: current directory)
- `--tolerance`: Convergence tolerance (default: 1e-5 for Fig4)
- `--max-iterations`: Maximum SCA iterations (default: 4000 for Fig4)

Run `python proposed_SCA.py --help` for complete options.

## Output Files

### Data Files
- `data_proposed_SCA_fig4.mat`: Complete results including:
  - SR_all: Sum rate values
  - CRB_all: CRB trace values  
  - Time_all: Computation times
  - Obj_all: Objective function values
  - Configuration parameters (antenna arrays, etc.)

### Plots
- `fig4_performance_vs_Nt.png`: Performance metrics vs number of antennas
- `fig4_fim_convergence.png`: Sensing Performance Convergence (Last Scenario)
- `fig4_objective_convergence.png`: Objective Function Convergence (Last Scenario)
- `fig4_sum_rate_convergence.png`: Sum Rate Convergence (Last Scenario)
- `plot_figure_OV_vs_Nt.png`: comparison with other algorithms, result of `plot_figure_OV_vs_Nt.py`

## Results Structure

### Main Arrays
All arrays have dimensions `(I_out, I_in)` where:
- I_out: Number of channel realizations (default: 100)
- I_in: Number of antenna configurations (5: for Nt = 8,16,32,64,128)

### Algorithm Details

1. **Antenna Scaling**: Tests Nt ∈ {8, 16, 32, 64, 128} antennas
2. **Channel Realizations**: Multiple random channel instances for each antenna configuration
3. **Objective Function**: `delta_c * sum_rate - delta_s * trace(inv(FIM))`
   - delta_c = 0.25 (fixed)
   - delta_s = 1 (fixed)

## Key Implementation Notes

### Antenna Array Scaling
- **Weight indexing**: weight = 1,2,3,4,5 corresponds to Nt = 8,16,32,64,128
- **Scaling formula**: Nt = Nth × Ntv = 4 × 2^(weight+1)
- **Receiver scaling**: Nr = Nrh × Nrv = 5 × 2^(weight+1)

### Beamforming Matrix Handling
- **Communication**: Wc = delta_c * H / ‖H‖ (per column normalization)
- **Sensing**: Random initialization for 6 sensing streams
- **Power constraint**: Normalized after each SCA iteration

### SCA Algorithm
- Same successive convex approximation as Figures 2 and 3
- Handles variable antenna array sizes in matrix operations
- Tighter convergence criteria (1e-5) and more iterations (4000)
- Uses same eigenvalue-based construction for dominant eigenvalue

## Expected Results

The analysis should show:
- Performance scaling with number of transmit antennas
- Computational complexity growth with antenna array size
- Comparison with other algorithms (WMMSE-SDR, FP-SGDA, LD Algorithm)
- Trade-offs between performance gains and computational cost

## Performance

Typical execution times (vary with antenna array size):
- Full run (100 realizations, all antenna sizes): ~30-60 minutes
- Quick test (10 realizations): ~3-6 minutes
- Per scenario: ~0.05-0.5 seconds (scales with antenna array size)

**⚠️ Performance Warning**: The current Python implementation has significant performance issues:
- **Actual runtime**: 18 hours, 7 minutes for full execution
- **Performance gap**: Extremely slow compared to original MATLAB version
- **Status**: Further investigation and optimization should be made

This performance degradation may be due to:
- Python vs MATLAB matrix operation efficiency differences
- Suboptimal NumPy/SciPy usage patterns
- Missing vectorization opportunities
- Memory allocation inefficiencies in tight loops

## Integration

This implementation integrates with:
- `utils/` module for shared functions (steering matrices, FIM calculation, etc.)
- `SimulationConfig` for parameter management
- Standard MATLAB .mat file format for compatibility
- Existing data files for algorithm comparison
- Plotting module `plot_figure_OV_vs_Nt.py` for visualization

## Algorithm Comparison

The Figure 4 analysis includes comparison with:
1. **WMMSE-SDR**: Weighted MMSE with semi-definite relaxation
2. **FP-SGDA**: Fractional programming with stochastic gradient descent
3. **Algorithm 1 (Ns=3M)**: Our SCA with full sensing streams
4. **Algorithm 1 (Ns=0)**: Our SCA with sensing-only optimization
5. **LD Algorithm 1**: Low-dimensional version of our SCA
