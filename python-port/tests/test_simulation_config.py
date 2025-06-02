#!/usr/bin/env python3
"""
Pytest tests for SimulationConfig class and demo_simulation_config.py

This module contains comprehensive tests for the simulation configuration
system, covering all functionality including default configs, custom configs,
scenario-specific configs, argument parsing, and helper methods.
"""

from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np
import argparse
import pytest
import sys
# Add the utils directory to the path
sys.path.append(str(Path(__file__).parent.parent / 'utils'))
# Import the modules to test
import demo_simulation_config
from simulation_config import SimulationConfig, db2pow

class TestSimulationConfig:
    """Test cases for SimulationConfig class"""

    @pytest.fixture
    def default_config(self):
        """Fixture providing a default configuration"""
        return SimulationConfig()

    def test_default_configuration(self, default_config):
        """Test that default configuration is created correctly"""
        config = default_config

        # Test antenna configuration
        assert config.Nth == 4
        assert config.Ntv == 4
        assert config.Nrh == 5
        assert config.Nrv == 4
        assert config.Nt == 16  # 4 * 4
        assert config.Nr == 20  # 5 * 4

        # Test system configuration
        assert config.K == 4
        assert config.K_max == 12
        assert config.M == 2
        assert config.M_max == 3

        # Test power configuration
        assert config.Pt_dBm == 10.0
        assert config.noise_c_dBm == 0.0
        assert config.noise_s_dBm == 0.0

        # Test computed power values
        assert abs(config.Pt - db2pow(10.0 - 30)) < 1e-10
        assert abs(config.noise_c - db2pow(0.0 - 30)) < 1e-10
        assert abs(config.noise_s - db2pow(0.0 - 30)) < 1e-10

        # Test sensing configuration
        assert config.L == 30
        assert abs(config.kappa - 2 * 30 / config.noise_s) < 1e-5

        # Test algorithm configuration
        assert config.tolerance == 1e-5
        assert config.max_iterations == 2000
        assert config.num_sensing_streams == config.Nt

        # Test weight parameters
        assert config.delta_s == 1.0
        assert config.delta_c == 0.25

    def test_custom_configuration(self):
        """Test custom configuration creation"""
        config = SimulationConfig(
            Nth=8, Ntv=8,
            K=6, M=3,
            Pt_dBm=15.0,
            L=64,
            delta_c=0.5
        )

        assert config.Nth == 8
        assert config.Ntv == 8
        assert config.Nt == 64  # 8 * 8
        assert config.K == 6
        assert config.M == 3
        assert config.Pt_dBm == 15.0
        assert config.L == 64
        assert config.delta_c == 0.5

        # Check computed values are updated
        assert abs(config.Pt - db2pow(15.0 - 30)) < 1e-10
        assert abs(config.kappa - 2 * 64 / config.noise_s) < 1e-5

    @pytest.mark.parametrize("db_val,expected", [
        (0, 1.0),
        (10, 10.0),
        (-10, 0.1),
        (30, 1000.0)
    ])
    def test_db2pow_function(self, db_val, expected):
        """Test db2pow utility function with various inputs"""
        assert abs(db2pow(db_val) - expected) < 1e-10

    @pytest.mark.parametrize("scenario,expected_params", [
        ('fig1', {'max_iterations': 2000, 'tolerance': 1e-5}),
        ('fig2', {'delta_all_step': 0.2, 'tolerance': 1e-4}),
        ('fig3', {'tolerance': 1e-4, 'L': 128, 'delta_c': 0.25}),
        ('fig4', {'tolerance': 1e-5, 'max_iterations': 4000}),
        ('fig5', {'tolerance': 1e-5, 'max_iterations': 4000}),
        ('fig6', {'delta_s': 0.0, 'tolerance': 1e-5, 'max_iterations': 4000}),
    ])
    def test_scenario_configurations(self, default_config, scenario, expected_params):
        """Test scenario-specific configurations"""
        scenario_config = default_config.get_scenario_config(scenario)

        for param, expected_value in expected_params.items():
            actual_value = getattr(scenario_config, param)
            assert actual_value == expected_value, f"Failed for {
                scenario}.{param}: {actual_value} != {expected_value}"

    def test_invalid_scenario_configuration(self, default_config):
        """Test invalid scenario returns unchanged config"""
        invalid_config = default_config.get_scenario_config('invalid_scenario')
        # Should return config with base parameters unchanged
        assert invalid_config.tolerance == default_config.tolerance

    def test_argument_parser_creation(self):
        """Test argument parser creation"""
        parser = SimulationConfig.create_argument_parser()
        assert isinstance(parser, argparse.ArgumentParser)

        # Test that parser has expected arguments
        help_text = parser.format_help()
        expected_args = ['--nth', '--ntv', '--k', '--m',
                         '--pt-dbm', '--delta-s', '--tolerance']
        for arg in expected_args:
            assert arg in help_text

    def test_from_args_with_mock_args(self):
        """Test configuration creation from parsed arguments"""
        # Create mock arguments
        mock_args = argparse.Namespace(
            nth=6, ntv=6, k=8, m=3,
            pt_dbm=12.0, delta_c=0.3,
            tolerance=1e-6, nrh=None, nrv=None
        )

        config = SimulationConfig.from_args(mock_args)

        assert config.Nth == 6
        assert config.Ntv == 6
        assert config.K == 8
        assert config.M == 3
        assert config.Pt_dBm == 12.0
        assert config.delta_c == 0.3
        assert config.tolerance == 1e-6

    def test_helper_methods(self):
        """Test helper methods"""
        config = SimulationConfig(random_seed=42)

        # Test generate_target_angles
        theta, phi = config.generate_target_angles()
        assert len(theta) == config.M
        assert len(phi) == config.M
        assert np.all(theta >= config.theta_range[0])
        assert np.all(theta <= config.theta_range[1])
        assert np.all(phi >= config.phi_range[0])
        assert np.all(phi <= config.phi_range[1])

        # Test reproducibility with same seed
        theta2, phi2 = config.generate_target_angles()
        np.testing.assert_array_almost_equal(theta, theta2)
        np.testing.assert_array_almost_equal(phi, phi2)

        # Test generate_alpha
        alpha = config.generate_alpha()
        assert len(alpha) == config.M
        assert np.all(np.iscomplex(alpha))

        # Test get_delta_range
        delta_range = config.get_delta_range()
        assert len(delta_range) > 0
        assert abs(delta_range[0] - config.delta_all_min) < 1e-5
        assert delta_range[-1] <= config.delta_all_max

    def test_string_representation(self, default_config):
        """Test __str__ method"""
        config_str = str(default_config)

        assert "Simulation Configuration:" in config_str
        assert f"Tx={default_config.Nt}" in config_str
        assert f"K={default_config.K}" in config_str
        assert f"Pt={default_config.Pt_dBm}dBm" in config_str

    def test_edge_cases(self):
        """Test edge cases and error conditions"""
        # Test with zero noise_s (should handle division by zero)
        config = SimulationConfig(noise_s_dBm=-float('inf'))
        # kappa should be 0 when noise_s is 0
        assert config.kappa == 0

        # Test with different M value for helper methods
        config = SimulationConfig(M=5, random_seed=123)
        theta, phi = config.generate_target_angles()
        alpha = config.generate_alpha()

        assert len(theta) == 5
        assert len(phi) == 5
        assert len(alpha) == 5


class TestDemoSimulationConfig:
    """Test cases for demo_simulation_config.py functions"""

    def test_demo_default_config(self, capsys):
        """Test demo_default_config function"""
        demo_simulation_config.demo_default_config()
        captured = capsys.readouterr()
        output = captured.out

        assert "=== Default Configuration ===" in output
        assert "Simulation Configuration:" in output
        assert "Antennas:" in output

    def test_demo_custom_config(self, capsys):
        """Test demo_custom_config function"""
        demo_simulation_config.demo_custom_config()
        captured = capsys.readouterr()
        output = captured.out

        assert "=== Custom Configuration ===" in output
        assert "Tx=64" in output  # Custom 8x8 antennas
        assert "K=6" in output    # Custom users

    def test_demo_scenario_configs(self, capsys):
        """Test demo_scenario_configs function"""
        demo_simulation_config.demo_scenario_configs()
        captured = capsys.readouterr()
        output = captured.out

        assert "=== Scenario-Specific Configurations ===" in output

        expected_scenarios = ['FIG1', 'FIG2', 'FIG3', 'FIG4', 'FIG5', 'FIG6']
        for scenario in expected_scenarios:
            assert f"{scenario} Configuration" in output

    @patch('demo_simulation_config.SimulationConfig.create_argument_parser')
    def test_demo_argument_parsing(self, mock_parser_creator, capsys):
        """Test demo_argument_parsing function"""
        # Create a mock parser and args
        mock_parser = MagicMock()
        mock_args = argparse.Namespace(
            nth=6, ntv=6, k=8, m=3, pt_dbm=12.0,
            delta_c=0.3, tolerance=1e-6
        )
        mock_parser.parse_args.return_value = mock_args
        mock_parser_creator.return_value = mock_parser

        demo_simulation_config.demo_argument_parsing()
        captured = capsys.readouterr()
        output = captured.out

        assert "=== Command Line Argument Parsing ===" in output
        assert "Configuration from command line arguments:" in output

    def test_demo_helper_methods(self, capsys):
        """Test demo_helper_methods function"""
        demo_simulation_config.demo_helper_methods()
        captured = capsys.readouterr()
        output = captured.out

        assert "=== Helper Methods ===" in output
        assert "Target angles" in output
        assert "Alpha values" in output
        assert "Delta range" in output

    @patch('demo_simulation_config.demo_default_config')
    @patch('demo_simulation_config.demo_custom_config')
    @patch('demo_simulation_config.demo_scenario_configs')
    @patch('demo_simulation_config.demo_argument_parsing')
    @patch('demo_simulation_config.demo_helper_methods')
    def test_main_function_success(self, mock_helper, mock_args,
                                   mock_scenario, mock_custom, mock_default):
        """Test main function successful execution"""
        result = demo_simulation_config.main()

        assert result == 0
        mock_default.assert_called_once()
        mock_custom.assert_called_once()
        mock_scenario.assert_called_once()
        mock_args.assert_called_once()
        mock_helper.assert_called_once()

    @patch('demo_simulation_config.demo_default_config')
    def test_main_function_with_exception(self, mock_default, capsys):
        """Test main function with exception handling"""
        mock_default.side_effect = Exception("Test error")

        result = demo_simulation_config.main()

        assert result == 1
        captured = capsys.readouterr()
        output = captured.out
        assert "Error during demonstration: Test error" in output


class TestIntegration:
    """Integration tests for the complete simulation config system"""

    def test_complete_workflow(self):
        """Test complete workflow from config creation to usage"""
        # Create configuration
        config = SimulationConfig(
            Nth=6, Ntv=6, K=8, M=3,
            Pt_dBm=12.0, L=64,
            random_seed=12345
        )

        # Test derived values
        assert config.Nt == 36
        assert abs(config.Pt - db2pow(12.0 - 30)) < 1e-10

        # Test scenario modification
        fig3_config = config.get_scenario_config('fig3')
        assert fig3_config.L == 128  # Should be overridden
        assert fig3_config.Nth == 6  # Should be preserved

        # Test helper methods work with modified config
        theta, phi = fig3_config.generate_target_angles()
        alpha = fig3_config.generate_alpha()
        delta_range = fig3_config.get_delta_range()

        assert len(theta) == 3  # M=3
        assert len(alpha) == 3  # M=3
        assert len(delta_range) > 0

    def test_command_line_integration(self):
        """Test command line argument integration"""
        test_args = [
            '--nth', '8', '--ntv', '4',
            '--k', '6', '--m', '3',
            '--pt-dbm', '15.0',
            '--tolerance', '1e-6',
            '--delta-c', '0.4'
        ]

        parser = SimulationConfig.create_argument_parser()
        args = parser.parse_args(test_args)
        config = SimulationConfig.from_args(args)

        assert config.Nth == 8
        assert config.Ntv == 4
        assert config.Nt == 32
        assert config.K == 6
        assert config.M == 3
        assert config.Pt_dBm == 15.0
        assert config.tolerance == 1e-6
        assert config.delta_c == 0.4

    def test_random_reproducibility(self):
        """Test that random number generation is reproducible"""
        config1 = SimulationConfig(random_seed=42, M=3)
        config2 = SimulationConfig(random_seed=42, M=3)

        theta1, phi1 = config1.generate_target_angles()
        alpha1 = config1.generate_alpha()

        theta2, phi2 = config2.generate_target_angles()
        alpha2 = config2.generate_alpha()

        np.testing.assert_array_almost_equal(theta1, theta2)
        np.testing.assert_array_almost_equal(phi1, phi2)
        np.testing.assert_array_almost_equal(alpha1, alpha2)


class TestErrorHandling:
    """Test error handling and edge cases"""

    def test_invalid_arguments(self):
        """Test handling of invalid arguments"""
        # Test with invalid scenario
        config = SimulationConfig()
        invalid_config = config.get_scenario_config('invalid_scenario')
        # Should return config without modification
        assert invalid_config.tolerance == config.tolerance

    def test_extreme_values(self):
        """Test with extreme parameter values"""
        # Test with very small/large values
        config = SimulationConfig(
            Nth=1, Ntv=1,  # Minimum antennas
            K=1, M=1,      # Minimum users/targets
            L=1,           # Minimum snapshots
            tolerance=1e-15  # Very small tolerance
        )

        assert config.Nt == 1
        assert config.K == 1
        assert config.M == 1
        assert config.L == 1

        # Test helper methods still work
        theta, phi = config.generate_target_angles()
        alpha = config.generate_alpha()

        assert len(theta) == 1
        assert len(alpha) == 1


# Pytest-specific test discovery and execution
if __name__ == "__main__":
    pytest.main([__file__])
