# Figure 5: Performance vs Number of Users (K)

This directory contains the Python translation of the MATLAB code for Figure 5, which analyzes the performance of the proposed SCA algorithm as a function of the number of communication users K.

## Files

- **`proposed_SCA.py`**: Main implementation of the SCA algorithm for Figure 5
- **`plot_figure.py`**: Script to generate the publication figure from saved data files
- **`data_*.mat`**: Pre-computed simulation results (MATLAB format)

## Algorithm Overview

The SCA (Successive Convex Approximation) algorithm optimizes the trade-off between communication sum rate and sensing performance (CRLB) by iteratively solving convex approximations of the original non-convex problem.

### Key Features:
- **Multi-user scenario**: Tests performance with K = 2, 4, 6, 8, 10, 12 communication users
- **ISAC system**: Joint optimization of communication and sensing beamforming
- **Convergence**: Uses tolerance-based stopping criteria for iterative optimization
- **Monte Carlo**: Averages results over multiple channel realizations

## Usage

### Full Simulation

Run the complete simulation with default parameters:

```bash
python proposed_SCA.py
```

### Custom Parameters

Run with custom simulation parameters:

```bash
python proposed_SCA.py --i-out 50 --i-in 6 --max-iterations 1000 --save-data --save-plots
```

### Generate Plots

Create publication-quality figures from existing data:

```bash
python plot_figure.py
```

## Parameters

### System Configuration
- **M = 2**: Number of sensing targets
- **K**: Number of communication users (varies: 2, 4, 6, 8, 10, 12)
- **Nt = 16**: Number of transmit antennas (4×4 array)
- **Nr = 20**: Number of receive antennas (5×4 array)

### Algorithm Configuration
- **δ_s = 1**: Sensing weight parameter
- **δ_c = 0.25**: Communication weight parameter
- **tolerance = 1e-5**: Convergence tolerance
- **max_iterations = 4000**: Maximum SCA iterations

### Simulation Configuration
- **I_out = 100**: Number of channel realizations
- **I_in = 6**: Number of K values to test
- **L = 30**: Number of sensing snapshots

## Output

The algorithm generates:

1. **Sum rate performance** for each K value
2. **CRLB trace** (sensing performance metric)
3. **Objective value** (weighted combination of sum rate and CRLB)
4. **Convergence data** for analysis
5. **Processing time** statistics

### Sample Output

```
Results summary (averaged over channels):
K       Objective       Sum Rate        CRLB Trace
--------------------------------------------------
2       -1.9613         5.9378          3.4457
4       -1.4258         9.2975          3.7502
6       -1.1496         11.1463         3.9361
8       -0.9535         12.4667         4.0701
10      -0.8274         13.2308         4.1352
12      -0.7114         14.0259         4.2179
```

## Command Line Options

Run `python proposed_SCA.py --help` to see all available options:

- **System parameters**: Antenna configurations, number of users/targets
- **Algorithm parameters**: Convergence criteria, iteration limits
- **Simulation parameters**: Monte Carlo iterations, random seeds
- **Output options**: Data saving, plot generation

## Data Files

The simulation saves results in MATLAB format compatible with the original implementation:

- **`data_SCA_Ns=3M.mat`**: Complete simulation results including:
  - `SR_all`: Sum rate results (I_out × I_in matrix)
  - `CRB_all`: CRLB trace results (I_out × I_in matrix)
  - `Time_all`: Processing time (I_out × I_in matrix)
  - `Obj_all`: Objective values (I_out × I_in matrix)
  - System parameters and configuration

## Performance Notes

- **Full simulation**: ~600 scenarios (100 channels × 6 K values)
- **Processing time**: ~0.4-0.5 seconds per scenario
- **Total runtime**: ~5-10 minutes for complete simulation
- **Memory usage**: Moderate (depends on antenna configuration)

## Validation

The Python implementation has been validated against the original MATLAB code:

1. **Convergence behavior**: Matches MATLAB implementation
2. **Final objective values**: Consistent results across implementations
3. **Processing time**: Comparable computational efficiency
4. **Output format**: Compatible data structures

## Troubleshooting

### Common Issues

1. **Import errors**: Ensure you're running from the correct directory and the virtual environment is activated
2. **Memory issues**: Reduce `I_out` or `I_in` for large-scale simulations
3. **Convergence problems**: Adjust `tolerance` or `max_iterations` parameters

### Performance Optimization

- Use `--i-out 10 --i-in 3` for quick testing
- Increase `--max-iterations` if convergence is slow
- Monitor memory usage for large antenna arrays
