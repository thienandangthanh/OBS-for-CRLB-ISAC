"""
Simulation Configuration Class for OBS-for-CRLB-ISAC

This module contains the SimulationConfig class that manages all environment
variables used in the MATLAB simulation files. It provides a centralized
configuration system with command-line argument parsing.

Extracted from MATLAB files:
- Fig1_convergence/propose_SCA.m
- Fig2_trade_off_region/proposed_SCA.m
- Fig3_performance_vs_Ns/proposed_SCA.m
- Fig_4performance_vs_Nt/proposed_SCA.m
- Fig_5performance_vs_K/proposed_SCA.m
- Fig_6performance_vs_Pt/proposed_SCA.m
"""

import argparse
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Tuple
import math
from utils.db2pow import db2pow


@dataclass
class SimulationConfig:
    """
    Configuration class containing all simulation environment variables.

    This class manages all parameters used across different MATLAB simulation
    scenarios and provides command-line argument parsing capabilities.
    """

    # === Antenna Configuration ===
    # Transmit antenna configuration
    Nth: int = 4  # Number of transmit antennas (horizontal)
    Ntv: int = 4  # Number of transmit antennas (vertical)

    # Receive antenna configuration
    Nrh: int = 5  # Number of receive antennas (horizontal)
    Nrv: int = 4  # Number of receive antennas (vertical)

    # Total antennas (computed automatically)
    Nt: int = field(init=False)  # Total transmit antennas
    Nr: int = field(init=False)  # Total receive antennas

    # === System Configuration ===
    K: int = 4          # Number of communication users/streams
    K_max: int = 12     # Maximum number of communication users
    M: int = 2          # Number of sensing targets
    M_max: int = 3      # Maximum number of targets

    # === Power and Noise Configuration ===
    Pt_dBm: float = 10.0      # Transmit power in dBm (before -30 offset)
    noise_c_dBm: float = 0.0  # Communication noise in dBm (before -30 offset)
    noise_s_dBm: float = 0.0  # Sensing noise in dBm (before -30 offset)

    # Power values in linear scale (computed automatically)
    Pt: float = field(init=False)      # Transmit power (linear)
    noise_c: float = field(init=False)  # Communication noise power (linear)
    noise_s: float = field(init=False)  # Sensing noise power (linear)

    # === Sensing Configuration ===
    # Number of sensing snapshots
    L: int = 30
    # Sensing parameter (computed as 2*L/noise_s)
    kappa: float = field(init=False)

    # === Algorithm Configuration ===
    tolerance: float = 1e-5         # Convergence tolerance
    max_iterations: int = 2000      # Maximum number of iterations
    num_sensing_streams: int = field(init=False)  # Number of sensing streams

    # === Weight/Trade-off Parameters ===
    delta_s: float = 1.0            # Sensing weight
    delta_c: float = 0.25           # Communication weight

    # === Monte Carlo Configuration ===
    channel_number: int = 50        # Number of channel realizations
    I_out: int = 50                 # Outer loop iterations
    I_in: int = 60                  # Inner loop parameter variations

    # === Range Parameters ===
    delta_all_min: float = -7.0     # Minimum delta value for trade-off
    delta_all_max: float = 4.8      # Maximum delta value for trade-off
    delta_all_step: float = 0.05    # Step size for delta sweep

    # === Random Seed ===
    random_seed: int = 0            # Random seed for reproducibility

    # === Angle Configuration ===
    theta_range: Tuple[float, float] = (-np.pi/3, np.pi/3)  # Theta angle range
    phi_range: Tuple[float, float] = (-np.pi/3, np.pi/3)    # Phi angle range

    # === Alpha (channel gain) Configuration ===
    alpha_base: float = 0.1         # Base value for alpha
    alpha_variance: float = 0.2     # Variance factor for alpha

    # === Variable Configuration for Different Figures ===
    # Fig 3: Ns performance (sensing streams variation)
    Ns_all: List[int] = field(default_factory=lambda: [
        0, 1, 2, 3, 4, 5, 6, 7, 8])

    # Fig 4: Nt performance (antenna variation)
    Nt_power_range: List[int] = field(
        default_factory=lambda: [2, 3, 4, 5, 6])  # 2^(i+1) for antenna count

    # Fig 5: K performance (user variation)
    K_multiplier_range: List[int] = field(
        default_factory=lambda: [1, 2, 3, 4, 5, 6])  # K = multiplier * 2

    # Fig 6: Pt performance (power variation)
    Pt_range_dBm: List[float] = field(
        default_factory=lambda: [-10, -5, 0, 5, 10, 15])  # Power range

    # === Convergence Figure Specific ===
    convergence_lin_values: List[float] = field(
        default_factory=lambda: [0.05, 0.1, 0.15])

    def __post_init__(self):
        """Compute derived values after initialization"""
        # Compute total antennas
        self.Nt = self.Nth * self.Ntv
        self.Nr = self.Nrh * self.Nrv

        # Compute power values (with -30dBm offset as in MATLAB)
        self.Pt = db2pow(self.Pt_dBm - 30)
        self.noise_c = db2pow(self.noise_c_dBm - 30)
        self.noise_s = db2pow(self.noise_s_dBm - 30)

        # Compute kappa
        self.kappa = 2 * self.L / self.noise_s if self.noise_s > 0 else 0

        # Set default number of sensing streams
        if not hasattr(self, 'num_sensing_streams') or self.num_sensing_streams is None:
            self.num_sensing_streams = self.Nt

    @classmethod
    def from_args(cls, args: Optional[argparse.Namespace] = None) -> 'SimulationConfig':
        """
        Create configuration from command line arguments.

        Args:
            args: Parsed arguments. If None, will parse from command line.

        Returns:
            SimulationConfig instance
        """
        if args is None:
            parser = cls.create_argument_parser()
            args = parser.parse_args()

        # Create config with parsed arguments
        config_dict = {}

        # Map argument names to config attributes
        arg_mapping = {
            'nth': 'Nth',
            'ntv': 'Ntv',
            'nrh': 'Nrh',
            'nrv': 'Nrv',
            'k': 'K',
            'k_max': 'K_max',
            'm': 'M',
            'm_max': 'M_max',
            'pt_dbm': 'Pt_dBm',
            'noise_c_dbm': 'noise_c_dBm',
            'noise_s_dbm': 'noise_s_dBm',
            'l': 'L',
            'tolerance': 'tolerance',
            'max_iterations': 'max_iterations',
            'delta_s': 'delta_s',
            'delta_c': 'delta_c',
            'channel_number': 'channel_number',
            'i_out': 'I_out',
            'i_in': 'I_in',
            'random_seed': 'random_seed',
            'alpha_base': 'alpha_base',
            'alpha_variance': 'alpha_variance'
        }

        # Extract arguments that are present
        for arg_name, config_attr in arg_mapping.items():
            if hasattr(args, arg_name) and getattr(args, arg_name) is not None:
                config_dict[config_attr] = getattr(args, arg_name)

        return cls(**config_dict)

    @staticmethod
    def create_argument_parser() -> argparse.ArgumentParser:
        """
        Create argument parser for simulation configuration.

        Returns:
            ArgumentParser instance
        """
        parser = argparse.ArgumentParser(
            description='OBS-for-CRLB-ISAC Simulation Configuration',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        # Antenna configuration
        antenna_group = parser.add_argument_group('Antenna Configuration')
        antenna_group.add_argument('--nth', type=int, default=4,
                                   help='Number of transmit antennas (horizontal)')
        antenna_group.add_argument('--ntv', type=int, default=4,
                                   help='Number of transmit antennas (vertical)')
        antenna_group.add_argument('--nrh', type=int, default=5,
                                   help='Number of receive antennas (horizontal)')
        antenna_group.add_argument('--nrv', type=int, default=4,
                                   help='Number of receive antennas (vertical)')

        # System configuration
        system_group = parser.add_argument_group('System Configuration')
        system_group.add_argument('--k', type=int, default=4,
                                  help='Number of communication users/streams')
        system_group.add_argument('--k-max', type=int, default=12,
                                  help='Maximum number of communication users')
        system_group.add_argument('--m', type=int, default=2,
                                  help='Number of sensing targets')
        system_group.add_argument('--m-max', type=int, default=3,
                                  help='Maximum number of targets')

        # Power and noise configuration
        power_group = parser.add_argument_group(
            'Power and Noise Configuration')
        power_group.add_argument('--pt-dbm', type=float, default=10.0,
                                 help='Transmit power in dBm')
        power_group.add_argument('--noise-c-dbm', type=float, default=0.0,
                                 help='Communication noise in dBm')
        power_group.add_argument('--noise-s-dbm', type=float, default=0.0,
                                 help='Sensing noise in dBm')

        # Sensing configuration
        sensing_group = parser.add_argument_group('Sensing Configuration')
        sensing_group.add_argument('--l', type=int, default=30,
                                   help='Number of sensing snapshots')

        # Algorithm configuration
        algo_group = parser.add_argument_group('Algorithm Configuration')
        algo_group.add_argument('--tolerance', type=float, default=1e-5,
                                help='Convergence tolerance')
        algo_group.add_argument('--max-iterations', type=int, default=2000,
                                help='Maximum number of iterations')

        # Weight parameters
        weight_group = parser.add_argument_group('Weight Parameters')
        weight_group.add_argument('--delta-s', type=float, default=1.0,
                                  help='Sensing weight')
        weight_group.add_argument('--delta-c', type=float, default=0.25,
                                  help='Communication weight')

        # Monte Carlo configuration
        mc_group = parser.add_argument_group('Monte Carlo Configuration')
        mc_group.add_argument('--channel-number', type=int, default=50,
                              help='Number of channel realizations')
        mc_group.add_argument('--i-out', type=int, default=50,
                              help='Outer loop iterations')
        mc_group.add_argument('--i-in', type=int, default=60,
                              help='Inner loop parameter variations')

        # Random seed
        parser.add_argument('--random-seed', type=int, default=0,
                            help='Random seed for reproducibility')

        # Alpha configuration
        alpha_group = parser.add_argument_group('Alpha Configuration')
        alpha_group.add_argument('--alpha-base', type=float, default=0.1,
                                 help='Base value for alpha')
        alpha_group.add_argument('--alpha-variance', type=float, default=0.2,
                                 help='Variance factor for alpha')

        return parser

    def get_delta_range(self) -> np.ndarray:
        """Get delta range for trade-off analysis"""
        return np.arange(self.delta_all_min, self.delta_all_max + self.delta_all_step,
                         self.delta_all_step)

    def generate_target_angles(self, rng: Optional[np.random.Generator] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate random target angles within specified ranges.

        Args:
            rng: Random number generator. If None, uses numpy default.

        Returns:
            Tuple of (theta, phi) arrays
        """
        if rng is None:
            rng = np.random.default_rng(self.random_seed)

        theta_min, theta_max = self.theta_range
        phi_min, phi_max = self.phi_range

        theta = theta_min + (theta_max - theta_min) * rng.random(self.M)
        phi = phi_min + (phi_max - phi_min) * rng.random(self.M)

        return theta, phi

    def generate_alpha(self, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """
        Generate random alpha (channel gain) values.

        Args:
            rng: Random number generator. If None, uses numpy default.

        Returns:
            Complex alpha array
        """
        if rng is None:
            rng = np.random.default_rng(self.random_seed)

        magnitude = self.alpha_base * \
            (1 + self.alpha_variance * rng.standard_normal(self.M))
        phase = 2 * np.pi * rng.random(self.M)

        return magnitude * np.exp(1j * phase)

    def get_scenario_config(self, scenario: str) -> 'SimulationConfig':
        """
        Get configuration for specific scenario (figure).

        Args:
            scenario: One of 'fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'fig6'

        Returns:
            Modified configuration for the scenario
        """
        # Create a copy of the current config by extracting only the init fields
        init_fields = {
            'Nth': self.Nth,
            'Ntv': self.Ntv,
            'Nrh': self.Nrh,
            'Nrv': self.Nrv,
            'K': self.K,
            'K_max': self.K_max,
            'M': self.M,
            'M_max': self.M_max,
            'Pt_dBm': self.Pt_dBm,
            'noise_c_dBm': self.noise_c_dBm,
            'noise_s_dBm': self.noise_s_dBm,
            'L': self.L,
            'tolerance': self.tolerance,
            'max_iterations': self.max_iterations,
            'delta_s': self.delta_s,
            'delta_c': self.delta_c,
            'channel_number': self.channel_number,
            'I_out': self.I_out,
            'I_in': self.I_in,
            'delta_all_min': self.delta_all_min,
            'delta_all_max': self.delta_all_max,
            'delta_all_step': self.delta_all_step,
            'random_seed': self.random_seed,
            'theta_range': self.theta_range,
            'phi_range': self.phi_range,
            'alpha_base': self.alpha_base,
            'alpha_variance': self.alpha_variance,
            'Ns_all': self.Ns_all,
            'Nt_power_range': self.Nt_power_range,
            'K_multiplier_range': self.K_multiplier_range,
            'Pt_range_dBm': self.Pt_range_dBm,
            'convergence_lin_values': self.convergence_lin_values
        }

        new_config = SimulationConfig(**init_fields)

        if scenario.lower() == 'fig1':
            # Convergence analysis
            new_config.max_iterations = 2000
            new_config.tolerance = 1e-5

        elif scenario.lower() == 'fig2':
            # Trade-off region
            new_config.delta_all_step = 0.2  # Coarser step
            new_config.tolerance = 1e-4

        elif scenario.lower() == 'fig3':
            # Performance vs Ns (sensing streams)
            new_config.tolerance = 1e-4
            new_config.L = 128  # Different L value
            new_config.delta_c = 0.25
            # Need to recompute dependent values
            new_config.__post_init__()

        elif scenario.lower() == 'fig4':
            # Performance vs Nt (antennas)
            new_config.tolerance = 1e-5
            new_config.max_iterations = 4000
            new_config.I_in = 5
            new_config.I_out = 100

        elif scenario.lower() == 'fig5':
            # Performance vs K (users)
            new_config.tolerance = 1e-5
            new_config.max_iterations = 4000

        elif scenario.lower() == 'fig6':
            # Performance vs Pt (power)
            new_config.delta_s = 0.0  # Communication only
            new_config.tolerance = 1e-5
            new_config.max_iterations = 4000

        return new_config

    def __str__(self) -> str:
        """String representation of configuration"""
        lines = [
            "Simulation Configuration:",
            f"  Antennas: Tx={self.Nt} ({self.Nth}x{self.Ntv}), Rx={
                self.Nr} ({self.Nrh}x{self.Nrv})",
            f"  System: K={self.K}, M={self.M}",
            f"  Power: Pt={self.Pt_dBm}dBm, Noise_c={
                self.noise_c_dBm}dBm, Noise_s={self.noise_s_dBm}dBm",
            f"  Sensing: L={self.L}, kappa={self.kappa:.2e}",
            f"  Weights: delta_s={self.delta_s}, delta_c={self.delta_c}",
            f"  Algorithm: tolerance={
                self.tolerance}, max_iter={self.max_iterations}",
            f"  Monte Carlo: channels={self.channel_number}, I_out={
                self.I_out}, I_in={self.I_in}"
        ]
        return "\n".join(lines)
