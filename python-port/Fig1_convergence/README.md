# Figure 1 Convergence Analysis

This directory contains the convergence analysis for the OBS-for-CRLB-ISAC system.
The code uses the centralized `SimulationConfig` class for environment variable management.

## Files

- `propose_SCA.py` - Main script for convergence analysis
- `plotfigure1_convergence.py` - Script to generate the combined convergence plot
- `example_usage.py` - Example usage demonstrating different configurations
- `README.md` - This documentation file

## Key Improvements

### 1. Centralized Configuration Management
- All environment variables are now managed through the `SimulationConfig` class
- No more hardcoded values scattered throughout the code
- Configuration can be easily modified through command-line arguments

### 2. Command-Line Interface
The script now supports comprehensive command-line arguments for all configuration parameters:

```bash
python propose_SCA.py --help
```

### 3. Enhanced Output Management
- Support for saving plots as files instead of displaying them
- Configurable output directory
- Better progress reporting and convergence information

### 4. Improved Code Structure
- Clear separation between configuration and algorithm logic
- Better error handling and user feedback
- More modular and maintainable code

## Usage

### Basic Usage

Run with default configuration:
```bash
python propose_SCA.py
```

### Advanced Usage

#### Save plots to files
```bash
python propose_SCA.py --save-plots --output-dir ./results
```

#### Modify system parameters
```bash
python propose_SCA.py --nth 6 --ntv 6 --k 6 --m 3
```

#### Change power settings
```bash
python propose_SCA.py --pt-dbm 15 --noise-c-dbm -5 --noise-s-dbm -5
```

#### Adjust algorithm parameters
```bash
python propose_SCA.py --tolerance 1e-6 --max-iterations 3000
```

#### Combine multiple options
```bash
python propose_SCA.py \
    --nth 8 --ntv 8 \
    --k 8 --m 3 \
    --pt-dbm 20 \
    --tolerance 1e-6 \
    --max-iterations 5000 \
    --save-plots \
    --output-dir ./high_performance_results
```

### Configuration Parameters

#### Antenna Configuration
- `--nth`: Number of transmit antennas (horizontal) [default: 4]
- `--ntv`: Number of transmit antennas (vertical) [default: 4]
- `--nrh`: Number of receive antennas (horizontal) [default: 5]
- `--nrv`: Number of receive antennas (vertical) [default: 4]

#### System Configuration
- `--k`: Number of communication users/streams [default: 4]
- `--m`: Number of sensing targets [default: 2]
- `--k-max`: Maximum number of communication users [default: 12]
- `--m-max`: Maximum number of targets [default: 3]

#### Power and Noise Configuration
- `--pt-dbm`: Transmit power in dBm [default: 10.0]
- `--noise-c-dbm`: Communication noise in dBm [default: 0.0]
- `--noise-s-dbm`: Sensing noise in dBm [default: 0.0]

#### Algorithm Configuration
- `--tolerance`: Convergence tolerance [default: 1e-5]
- `--max-iterations`: Maximum number of iterations [default: 2000]
- `--l`: Number of sensing snapshots [default: 30]

#### Weight Parameters
- `--delta-s`: Sensing weight [default: 1.0]
- `--delta-c`: Communication weight [default: 0.25]

#### Other Parameters
- `--random-seed`: Random seed for reproducibility [default: 0]
- `--alpha-base`: Base value for alpha [default: 0.1]
- `--alpha-variance`: Variance factor for alpha [default: 0.2]

### Output Files

The script generates the following output files:

1. **data_convergence.mat** - MATLAB format file containing convergence data
2. **fig1_sum_rate_vs_iteration.png** - Sum Rate vs Iteration plot (if --save-plots is used)
3. **fig1_inverse_fim_trace_vs_iteration.png** - Trace of inverse FIM vs Iteration plot (if --save-plots is used)
4. **fig1_objective_value_vs_iteration.png** - Total objective value vs Iteration plot (if --save-plots is used)

### Additional Plotting

Use the `plotfigure1_convergence.py` script to generate the main publication-quality convergence plot:

```bash
# Generate and display the combined convergence plot
python plotfigure1_convergence.py

# Save the plot without displaying it
python plotfigure1_convergence.py --save-only --output-dir ./results

# Custom output file name and directory
python plotfigure1_convergence.py --output-file custom_convergence.png --output-dir ./figures --save-only
```

This script generates:
- **fig1_convergence_combined.png** - Combined convergence plot showing all three δc values on a single figure

## Examples

Run the example script to see different usage scenarios:

```bash
python example_usage.py
```

This will run several examples with different configurations and save results to separate directories.
