# Figure 3: Performance vs Number of Sensing Streams (Ns)

This directory contains the Python translation of the MATLAB Figure 3 analysis for the OBS-for-CRLB-ISAC system.

## Overview

Figure 3 analyzes how the system objective function varies with the number of sensing streams (Ns) for different numbers of communication users (K).
Unlike Figure 2 which analyzes trade-offs between communication and sensing weights, Figure 3 focuses on the impact of dedicating different numbers of streams to sensing.

## Implementation Scope

This implementation currently includes:
- The main SCA algorithm for K=2,3,4 users (proposed_SCA.py)
- The sensing-only case (K=0) will be implemented separately from proposed__SCA_sensing_only.m
- Data and result from WMMSE-SDR algorithm are preserved as-is from MATLAB implementation
- Only the SCA algorithm is currently ported to Python

## Key Differences from Figure 2 SCA Implementation

The following differences are specific to the SCA algorithm implementations between Figure 2 and Figure 3:

| Parameter | Figure 2 SCA | Figure 3 SCA |
|-----------|----------|----------|
| **M_max** | 2 | 3 |
| **m_target** | 2 | 3 |
| **L** | 30 | 128 |
| **num_sensing_streams** | Nt (fixed) | weight-1 (varies: 0-8) |
| **delta_c** | 10^(delta_all(weight)) | 0.25 (fixed) |
| **delta_s** | 1 (fixed) | 1 (fixed) |
| **User loop** | None | 2, 3, 4 users |
| **H_all dimensions** | (I_out, Nt, K) | (Nt, K_max, I_out) |
| **Wc initialization** | Random | delta_c*H/||H|| (per column normalization) |
| **Main output** | SR_all, CRB_all, Time_all | Obj_all |

## Usage

### Basic Usage
```bash
python proposed_SCA.py
```

### Quick Test (Reduced Parameters)
```bash
python proposed_SCA.py --i-out 5 --save-data --save-plots
```

### Custom Configuration
```bash
python proposed_SCA.py --l 64 --tolerance 1e-5 --max-iterations 1000 --output-dir ./results
```

## Command Line Options

All standard simulation parameters are available via command line. Key options:

- `--i-out`: Number of channel realizations (default: 50)
- `--save-data`: Save results to data files
- `--save-plots`: Save plots as files instead of displaying
- `--output-dir`: Directory for output files (default: current directory)
- `--l`: Number of sensing snapshots (default: 128 for Fig3)
- `--tolerance`: Convergence tolerance (default: 1e-4 for Fig3)

Run `python proposed_SCA.py --help` for complete options.

## Output Files

### Data Files
- `data_proposed_SCA_fig3.mat`: Complete results including:
  - Obj_all: Main objective values
  - SR_all: Sum rate values
  - CRB_all: CRB trace values
  - Time_all: Computation times
  - Configuration parameters (Nth, Ntv, Nt, etc.)
- `data.mat`: Intermediate results (Obj_all only, saved after each user)

### Plots
- `fig3_performance_vs_Ns.png`: Objective value vs number of sensing streams

## Results Structure

### Obj_all
3D array with dimensions `(4, I_out, I_in)`:
- Dimension 1: User configurations (0-3 indices for 1-4 users)
  - Index 0: Reserved for K=0 case (to be implemented from proposed__SCA_sensing_only.m)
  - Indices 1-3: For K=1,2,3 communication users (2,3,4 total users)
- Dimension 2: Channel realizations (I_out)
- Dimension 3: Number of sensing streams (0-8)

### Algorithm Details

1. **User Loop**: Tests K=1, K=2, K=3 communication users (2,3,4 total users)
   - K=0 case will be implemented separately from proposed__SCA_sensing_only.m
2. **Sensing Streams**: Varies from 0 to 8 streams (Ns_all=[0,1,2,3,4,5,6,7,8])
3. **Channel Realizations**: Multiple random channel instances
4. **Objective Function**: `delta_c * sum_rate - delta_s * trace(inv(FIM))`
   - delta_c = 0.25 (fixed)
   - delta_s = 1 (fixed)

## Key Implementation Notes

### Sensing Stream Handling
- When `num_sensing_streams = 0`: Only communication beamforming (Wc)
- When `num_sensing_streams > 0`: Combined communication + sensing beamforming

### Beamforming Initialization
- **Communication**: `Wc = delta_c * H / ||H||` (per column normalization)
- **Sensing**: Random initialization when `num_sensing_streams > 0`

### SCA Algorithm
- Same successive convex approximation as Figure 2
- Handles variable number of sensing streams in matrix construction
- Power constraint normalization after each iteration
- Convergence based on norm of beamforming matrix difference

## Expected Results

The plot should show:
- Objective function values vs number of sensing streams (0-8)
- Separate curves for K=1, K=2, K=3 users (2, 3, 4 total users)
- Generally increasing performance with more sensing streams
- Trade-offs between communication and sensing requirements

## Performance

Typical execution times:
- Full run (50 realizations): ~10-15 minutes
- Quick test (5 realizations): ~1-2 minutes
- Per scenario: ~0.1-0.2 seconds

## Integration

This implementation integrates with:
- `utils/` module for shared functions
- `SimulationConfig` for parameter management
- Standard MATLAB .mat file format for compatibility
- Optional plotting module `plot_SCA_performance_vs_Ns.py`

## Future Work

1. Implementation of sensing-only case (K=0) from proposed__SCA_sensing_only.m
2. Extract core SCA algorithm for better code maintenance
