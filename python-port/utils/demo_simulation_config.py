#!/usr/bin/env python3
"""
Demonstration script for SimulationConfig class

This script shows how to use the SimulationConfig class for managing
simulation parameters in the OBS-for-CRLB-ISAC project.
"""

import sys
import argparse
from simulation_config import SimulationConfig


def demo_default_config():
    """Demonstrate default configuration"""
    print("=== Default Configuration ===")
    config = SimulationConfig()
    print(config)
    print()


def demo_custom_config():
    """Demonstrate custom configuration"""
    print("=== Custom Configuration ===")
    config = SimulationConfig(
        Nth=8, Ntv=8,  # Larger antenna array
        K=6, M=3,      # More users and targets
        Pt_dBm=15.0,   # Higher power
        L=64,          # More snapshots
        delta_c=0.5    # Different weight
    )
    print(config)
    print()


def demo_scenario_configs():
    """Demonstrate scenario-specific configurations"""
    print("=== Scenario-Specific Configurations ===")

    base_config = SimulationConfig()

    scenarios = ['fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'fig6']

    for scenario in scenarios:
        print(f"\n--- {scenario.upper()} Configuration ---")
        scenario_config = base_config.get_scenario_config(scenario)
        print(f"Tolerance: {scenario_config.tolerance}")
        print(f"Max iterations: {scenario_config.max_iterations}")
        print(f"Delta_s: {scenario_config.delta_s}")
        print(f"Delta_c: {scenario_config.delta_c}")
        if hasattr(scenario_config, 'L'):
            print(f"L: {scenario_config.L}")


def demo_argument_parsing():
    """Demonstrate command line argument parsing"""
    print("=== Command Line Argument Parsing ===")

    # Simulate command line arguments
    test_args = [
        '--nth', '6',
        '--ntv', '6',
        '--k', '8',
        '--m', '3',
        '--pt-dbm', '12.0',
        '--delta-c', '0.3',
        '--tolerance', '1e-6'
    ]

    # Parse arguments
    parser = SimulationConfig.create_argument_parser()
    args = parser.parse_args(test_args)

    # Create config from arguments
    config = SimulationConfig.from_args(args)
    print("Configuration from command line arguments:")
    print(config)
    print()


def demo_helper_methods():
    """Demonstrate helper methods"""
    print("=== Helper Methods ===")

    config = SimulationConfig(random_seed=42)

    # Generate target angles
    theta, phi = config.generate_target_angles()
    print(f"Target angles - Theta: {theta}, Phi: {phi}")

    # Generate alpha values
    alpha = config.generate_alpha()
    print(f"Alpha values: {alpha}")

    # Get delta range
    delta_range = config.get_delta_range()
    print(f"Delta range: {len(delta_range)} points from {
          delta_range[0]:.2f} to {delta_range[-1]:.2f}")
    print()


def main():
    """Main demonstration function"""
    print("SimulationConfig Class Demonstration")
    print("=" * 50)
    print()

    try:
        demo_default_config()
        demo_custom_config()
        demo_scenario_configs()
        demo_argument_parsing()
        demo_helper_methods()

        print("Demo completed successfully!")

    except Exception as e:
        print(f"Error during demonstration: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
