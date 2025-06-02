# SimulationConfig Class Documentation

## Overview

The `SimulationConfig` class provides a centralized configuration system for the OBS-for-CRLB-ISAC simulation project. It extracts and manages all environment variables used across different MATLAB simulation files, providing a unified Python interface with command-line argument parsing capabilities.

## Key Features

- **Centralized Configuration**: All simulation parameters in one class
- **Command-line Interface**: Full argument parsing support
- **Scenario-specific Configurations**: Pre-configured setups for different figures
- **Automatic Derived Values**: Computed fields like total antennas and power conversions
- **Helper Methods**: Generate random angles, alpha values, and parameter ranges

## Environment Variables Extracted

Based on analysis of the MATLAB files, the following environment variables were extracted:

### Antenna Configuration
- `Nth`, `Ntv`: Transmit antenna dimensions (horizontal × vertical)
- `Nrh`, `Nrv`: Receive antenna dimensions (horizontal × vertical)
- `Nt`, `Nr`: Total transmit/receive antennas (computed)

### System Configuration
- `K`: Number of communication users/streams
- `K_max`: Maximum number of communication users  
- `M`: Number of sensing targets
- `M_max`: Maximum number of targets

### Power and Noise Configuration
- `Pt_dBm`: Transmit power in dBm (with -30dBm offset applied)
- `noise_c_dBm`: Communication noise power in dBm
- `noise_s_dBm`: Sensing noise power in dBm
- Linear power values are computed automatically

### Algorithm Configuration
- `tolerance`: Convergence tolerance
- `max_iterations`: Maximum number of iterations
- `L`: Number of sensing snapshots
- `kappa`: Sensing parameter (computed as 2*L/noise_s)

### Weight Parameters
- `delta_s`: Sensing weight
- `delta_c`: Communication weight

### Monte Carlo Configuration
- `channel_number`: Number of channel realizations
- `I_out`: Outer loop iterations
- `I_in`: Inner loop parameter variations

## Differences Between MATLAB Files

After analyzing all the propose_SCA.m/proposed_SCA.m files across different figures, the main differences found are:

### Figure 1 (Convergence Analysis)
- Standard parameters with `tolerance=1e-5`
- `max_iterations=2000`

### Figure 2 (Trade-off Region)
- Coarser delta step size: `delta_all_step=0.2`
- Relaxed tolerance: `tolerance=1e-4`

### Figure 3 (Performance vs Ns)
- Different L value: `L=128` vs standard `L=30`
- `tolerance=1e-4`
- Focus on sensing streams variation

### Figure 4 (Performance vs Nt)
- Higher iteration limit: `max_iterations=4000`
- `tolerance=1e-5`
- Antenna count variations

### Figure 5 (Performance vs K)  
- Higher iteration limit: `max_iterations=4000`
- `tolerance=1e-5`
- User count variations: `K = weight*2`

### Figure 6 (Performance vs Pt)
- Communication-only mode: `delta_s=0.0`
- Power variations: `Pt=db2pow(5*weight-10-30)`
- `max_iterations=4000`

### Common Variables Across All Files
- `Nth=4, Ntv=4` (transmit antennas)
- `Nrh=5, Nrv=4` (receive antennas)
- `M=2` (sensing targets)
- `Pt=db2pow(10-30)` (base power level)
- `noise_c=db2pow(0-30), noise_s=db2pow(0-30)` (noise levels)
- `alpha=0.1*(1+0.2*randn(M,1)).*exp(1j*2*pi*rand(M,1))` (channel gains)
- `theta, phi` angles in range `[-π/3, π/3]`

## Usage Examples

### Basic Usage

```python
from utils import SimulationConfig

# Create default configuration
config = SimulationConfig()
print(config)

# Create custom configuration  
config = SimulationConfig(
    Nth=8, Ntv=8,    # Larger antenna array
    K=6, M=3,        # More users and targets
    Pt_dBm=15.0,     # Higher power
    delta_c=0.5      # Different weight
)
```

### Command Line Usage

```bash
# Run with custom parameters
python my_simulation.py --nth 8 --ntv 8 --k 6 --m 3 --pt-dbm 15.0 --delta-c 0.5

# See all available options
python my_simulation.py --help
```

### Scenario-specific Configurations

```python
# Get configuration for specific figure
base_config = SimulationConfig()
fig3_config = base_config.get_scenario_config('fig3')
print(f"Fig 3 uses L={fig3_config.L} snapshots")
```

### Helper Methods

```python
config = SimulationConfig(random_seed=42)

# Generate target angles
theta, phi = config.generate_target_angles()

# Generate alpha channel gains
alpha = config.generate_alpha()

# Get delta range for trade-off analysis
delta_range = config.get_delta_range()
```

### Integration in Your Code

```python
import argparse
from utils import SimulationConfig

def main():
    # Parse command line arguments
    config = SimulationConfig.from_args()
    
    # Use configuration in your simulation
    print(f"Running simulation with {config.Nt} transmit antennas")
    print(f"Power: {config.Pt:.2e}, Noise: {config.noise_s:.2e}")
    
    # Generate simulation parameters
    theta, phi = config.generate_target_angles()
    alpha = config.generate_alpha()
    
    # Your simulation code here...

if __name__ == "__main__":
    main()
```

## File Structure

```
python-port/utils/
├── simulation_config.py      # Main configuration class
├── demo_simulation_config.py # Demonstration script
└── __init__.py              # Module exports
```

## Dependencies

- `numpy`: For numerical operations and random number generation
- `argparse`: For command-line argument parsing (standard library)
- `dataclasses`: For structured data (standard library)
- `typing`: For type hints (standard library) 