# Python Port of ISAC Beamforming Simulation Code

This directory contains Python translations of the MATLAB simulation code for the paper "Optimal ISAC Beamforming Structure and Efficient Algorithms for Sum Rate and CRLB Balancing". The original MATLAB code can be found in the parent directory `OBS_for_CRLB_ISAC/`.

## About This Project

This Python port aims to provide equivalent functionality to the original MATLAB simulations for researchers who prefer working in Python or need to integrate the algorithms into Python-based workflows. The code generates publication-quality figures demonstrating:

- **Figure 1**: Convergence behavior of optimization algorithms
- **Figure 2**: Trade-off regions between communication and sensing performance
- **Figure 3**: Performance vs. number of sensing streams
- **Figure 4**: Performance vs. number of transmit antennas
- **Figure 5**: Performance vs. number of users
- **Figure 6**: Performance vs. transmit power

## Project Structure

```
python-port/
├── README.md
├── requirements.txt
├── Fig1_convergence/
├── Fig2_trade_off_region/
├── Fig3_performance_vs_Ns/
└── ... (other figures)
```

## Requirements

Before running this project, ensure you have:

- **Python 3.7 or higher** installed on your system
- **pip** package manager
- The corresponding `.mat` data files (copied from the original MATLAB directories)

## Setting Up the Virtual Environment

### 1. Create a Virtual Environment

Navigate to the `python-port` directory and create a virtual environment:

```bash
# Create virtual environment named 'venv'
python3 -m venv venv
```

### 2. Activate the Virtual Environment

**On Linux/macOS:**
```bash
source venv/bin/activate
```
With fish shell, use this instead
```fish
source venv/bin/activate.fish
```

**On Windows:**
```bash
venv\Scripts\activate
```

You should see `(venv)` prefix in your terminal prompt indicating the virtual environment is active.

### 3. Install Dependencies

Install the required Python packages using pip:

```bash
pip install -r requirements.txt
```

## Running Individual Figure Scripts

Each figure has its own subdirectory with the corresponding Python script:

```bash
# Example: Run Figure 1 convergence plot
cd Fig1_convergence
python plotfigure1_convergence.py

# Example: Run Figure 2 trade-off region
cd Fig2_trade_off_region
python plotfigure2_tradeoff_region.py
```

## Data Files

The Python scripts use the same `.mat` data files as the original MATLAB code. These files contain pre-computed simulation results and should be copied from the corresponding MATLAB directories:

- `Fig1_convergence/data_convergence.mat`
- `Fig2_trade_off_region/data_*.mat`
- `Fig3_performance_vs_Ns/data_*.mat`
- And so on...

## Translation Status

- ✅ **Figure 1 (Convergence)**: Complete
- ✅ **Figure 2 (Trade-off Region)**: Complete
- ✅ **Figure 3 (Performance vs Ns)**: Complete
- ✅ **Figure 4 (Performance vs Nt)**: Complete
- ✅ **Figure 5 (Performance vs K)**: Complete
- ✅ **Figure 6 (Performance vs Pt)**: Complete

## Deactivating the Virtual Environment

When you're done working on the project, deactivate the virtual environment:

```bash
deactivate
```

## Original Paper Reference

This Python port is based on the MATLAB code for:

**"Optimal ISAC Beamforming Structure and Efficient Algorithms for Sum Rate and CRLB Balancing"**

ArXiv: https://arxiv.org/pdf/2503.09489

## Troubleshooting

### Common Issues:

1. **ModuleNotFoundError**: Make sure the virtual environment is activated and dependencies are installed
2. **File not found error**: Ensure the corresponding `.mat` files are copied to the correct subdirectories
3. **Data structure issues**: The scripts include debugging output to help identify MATLAB data structure problems

### Getting Help:

- Check that you're running scripts from the correct subdirectory
- Ensure all required `.mat` files are present
- Verify that the virtual environment is activated before running scripts
- Use the debugging output in scripts to inspect data structures if needed
