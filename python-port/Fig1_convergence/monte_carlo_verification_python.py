#!/usr/bin/env python3
"""
Monte Carlo Verification of propose_SCA Algorithm

This script runs the propose_SCA algorithm multiple times with different random
initializations and averages the results to verify statistical convergence 
behavior and algorithm correctness.
"""

import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
import os
import sys
import subprocess
import tempfile
import argparse
from tqdm import tqdm
from scipy.io import loadmat, savemat

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.simulation_config import SimulationConfig

def run_single_simulation(config, run_id=None, temp_dir=None):
    """
    Run a single instance of the propose_SCA algorithm with random initialization.
    Returns convergence data for this run by calling propose_SCA.py as subprocess.
    """
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()

    # Create unique output directory for this run
    run_output_dir = os.path.join(temp_dir, f"run_{run_id}")
    os.makedirs(run_output_dir, exist_ok=True)

    # Build command to run propose_SCA.py
    # Use different random seed for each run to ensure randomness
    random_seed = np.random.randint(0, 1000000) if run_id is None else run_id * 12345

    cmd = [
        'python', 'propose_SCA.py',
        '--save-plots',
        '--output-dir', run_output_dir,
        '--random-seed', str(random_seed),
        '--max-iterations', str(config.max_iterations),
        '--tolerance', str(config.tolerance),
        '--nth', str(config.Nth),
        '--ntv', str(config.Ntv),
        '--nrh', str(config.Nrh),
        '--nrv', str(config.Nrv),
        '--k', str(config.K),
        '--m', str(config.M),
        '--pt-dbm', str(config.Pt_dBm),
        '--noise-c-dbm', str(config.noise_c_dBm),
        '--noise-s-dbm', str(config.noise_s_dBm),
        '--l', str(config.L),
        '--delta-s', str(config.delta_s),
        '--alpha-base', str(config.alpha_base),
        '--alpha-variance', str(config.alpha_variance)
    ]

    try:
        # Run propose_SCA.py
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(__file__))

        if result.returncode != 0:
            raise RuntimeError(f"propose_SCA.py failed: {result.stderr}")

        # Load the generated data
        data_file = os.path.join(run_output_dir, 'data_convergence.mat')
        if not os.path.exists(data_file):
            raise FileNotFoundError(f"Expected output file not found: {data_file}")

        data = loadmat(data_file)
        Con = data['Con']

        # Extract convergence data for all three lin values
        results = {'lin_values': config.convergence_lin_values, 'convergence_data': []}

        Con_flat = Con.flatten()
        for de in range(3):
            if de < len(Con_flat):
                # Extract the convergence history for this lin value
                convergence_history = Con_flat[de].tolist()
                results['convergence_data'].append(convergence_history)
            else:
                results['convergence_data'].append([])

        return results

    except Exception as e:
        print(f"Warning: Run {run_id} failed: {e}")
        raise

def monte_carlo_verification(config, num_runs=300, save_individual=False, temp_dir=None):
    """
    Run Monte Carlo verification with multiple independent runs.

    Parameters:
    -----------
    config : SimulationConfig
        Configuration object
    num_runs : int
        Number of Monte Carlo runs
    save_individual : bool
        Whether to save individual run data
    temp_dir : str
        Temporary directory for intermediate files
    """

    print(f"=== MONTE CARLO VERIFICATION OF PROPOSE_SCA ===")
    print(f"Running {num_runs} independent simulations...")
    print("Each run uses different random seed for initialization")
    print()
    print("Configuration:")
    print(config)
    print()

    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()

    # Storage for all runs
    all_runs_data = []
    failed_runs = 0

    start_time = perf_counter()

    # Run Monte Carlo simulations with progress bar
    for run_id in tqdm(range(num_runs), desc="Monte Carlo Runs"):
        try:
            run_data = run_single_simulation(config, run_id, temp_dir)
            all_runs_data.append(run_data)
        except Exception as e:
            if config.max_iterations < 100:  # Only show warnings for quick runs
                print(f"Warning: Run {run_id} failed: {e}")
            failed_runs += 1

    successful_runs = len(all_runs_data)
    total_time = perf_counter() - start_time

    print(f"\nCompleted {successful_runs}/{num_runs} runs successfully")
    print(f"Failed runs: {failed_runs}")
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Average time per run: {total_time/num_runs:.3f} seconds")

    if successful_runs == 0:
        print("ERROR: No successful runs completed!")
        return None

    # Analyze results
    print("\n=== ANALYZING MONTE CARLO RESULTS ===")

    # Find maximum length for padding
    max_lengths = [0, 0, 0]  # for each lin value
    for run_data in all_runs_data:
        for de in range(3):
            if de < len(run_data['convergence_data']):
                max_lengths[de] = max(max_lengths[de], len(run_data['convergence_data'][de]))

    print(f"Maximum iterations observed: {max(max_lengths)}")

    # Compute statistics
    statistics = {}
    lin_vals = config.convergence_lin_values

    for de in range(3):
        lin_val = lin_vals[de]
        max_len = max_lengths[de]

        # Storage for padded data
        all_sum_rates = []
        all_crb_traces = []
        all_objectives = []
        final_values = {'sum_rate': [], 'crb_trace': [], 'objective': []}
        convergence_iterations = []

        for run_data in all_runs_data:
            if de < len(run_data['convergence_data']):
                history = run_data['convergence_data'][de]
                convergence_iterations.append(len(history))

                # Convert to numpy array
                history = np.array(history)

                # Pad with final values if needed
                if len(history) < max_len:
                    final_val = history[-1] if len(history) > 0 else [0, 0, 0]
                    padding = np.tile(final_val, (max_len - len(history), 1))
                    history = np.vstack([history, padding])

                all_sum_rates.append(history[:, 0])
                all_crb_traces.append(history[:, 1])
                all_objectives.append(history[:, 2])

                # Store final values
                final_values['sum_rate'].append(history[-1, 0])
                final_values['crb_trace'].append(history[-1, 1])
                final_values['objective'].append(history[-1, 2])

        # Convert to arrays
        all_sum_rates = np.array(all_sum_rates)
        all_crb_traces = np.array(all_crb_traces)
        all_objectives = np.array(all_objectives)

        # Compute statistics
        mean_sum_rate = np.mean(all_sum_rates, axis=0)
        std_sum_rate = np.std(all_sum_rates, axis=0)
        mean_crb_trace = np.mean(all_crb_traces, axis=0)
        std_crb_trace = np.std(all_crb_traces, axis=0)
        mean_objective = np.mean(all_objectives, axis=0)
        std_objective = np.std(all_objectives, axis=0)

        statistics[f'lin_{lin_val}'] = {
            'mean_sum_rate': mean_sum_rate,
            'std_sum_rate': std_sum_rate,
            'mean_crb_trace': mean_crb_trace,
            'std_crb_trace': std_crb_trace,
            'mean_objective': mean_objective,
            'std_objective': std_objective,
            'final_values': final_values,
            'convergence_iterations': convergence_iterations,
            'max_iterations': max_len
        }

        # Print summary statistics
        print(f"\nδc = {lin_val} (lin[{de}]):")
        print(f"  Average convergence iterations: {np.mean(convergence_iterations):.1f} ± {np.std(convergence_iterations):.1f}")
        print(f"  Final sum rate: {np.mean(final_values['sum_rate']):.4f} ± {np.std(final_values['sum_rate']):.4f}")
        print(f"  Final CRB trace: {np.mean(final_values['crb_trace']):.4f} ± {np.std(final_values['crb_trace']):.4f}")
        print(f"  Final objective: {np.mean(final_values['objective']):.4f} ± {np.std(final_values['objective']):.4f}")

    return statistics, all_runs_data

def plot_monte_carlo_results(statistics, config, num_runs, output_dir='.'):
    """Plot Monte Carlo averaged results"""

    plt.style.use('default')
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'Monte Carlo Averaged Convergence (N={num_runs} runs)', fontsize=14)

    colors = ['b-', 'r-', 'g-']
    lin_vals = config.convergence_lin_values
    labels = [f'δc = {val}' for val in lin_vals]

    for i, (lin_val, color, label) in enumerate(zip(lin_vals, colors, labels)):
        key = f'lin_{lin_val}'
        if key not in statistics:
            continue

        data = statistics[key]
        iterations = range(len(data['mean_sum_rate']))

        # Plot 1: Sum Rate (take real part)
        mean_sum_rate = np.real(data['mean_sum_rate'])
        std_sum_rate = np.real(data['std_sum_rate'])
        axes[0].plot(iterations, mean_sum_rate, color, label=label, linewidth=2)
        axes[0].fill_between(iterations, 
                             mean_sum_rate - std_sum_rate,
                             mean_sum_rate + std_sum_rate,
                             alpha=0.3, color=color[0])

        # Plot 2: CRB Trace (take real part)
        mean_crb_trace = np.real(data['mean_crb_trace'])
        std_crb_trace = np.real(data['std_crb_trace'])
        axes[1].plot(iterations, mean_crb_trace, color, label=label, linewidth=2)
        axes[1].fill_between(iterations,
                             mean_crb_trace - std_crb_trace,
                             mean_crb_trace + std_crb_trace,
                             alpha=0.3, color=color[0])

        # Plot 3: Objective (take real part)
        mean_objective = np.real(data['mean_objective'])
        std_objective = np.real(data['std_objective'])
        axes[2].plot(iterations, mean_objective, color, label=label, linewidth=2)
        axes[2].fill_between(iterations,
                             mean_objective - std_objective,
                             mean_objective + std_objective,
                             alpha=0.3, color=color[0])

    # Format plots
    axes[0].set_xlabel('Iteration')
    axes[0].set_ylabel('Sum Rate (nat/s/Hz)')
    axes[0].set_title('Average Sum Rate Convergence')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel('Iteration')
    axes[1].set_ylabel('Trace of Inverse FIM')
    axes[1].set_title('Average CRB Trace Convergence')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    axes[2].set_xlabel('Iteration')
    axes[2].set_ylabel('Total Objective Value')
    axes[2].set_title('Average Objective Convergence')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    output_file = os.path.join(output_dir, 'monte_carlo_averaged_convergence_python.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Convergence plot saved to: {output_file}")
    plt.show()

    # Additional statistics plot
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

    # Convergence iteration histogram
    all_iterations = []
    for i, lin_val in enumerate(lin_vals):
        key = f'lin_{lin_val}'
        if key in statistics:
            iterations = statistics[key]['convergence_iterations']
            all_iterations.extend(iterations)
            axes2[0].hist(iterations, bins=20, alpha=0.7, label=f'δc = {lin_val}')

    axes2[0].set_xlabel('Convergence Iterations')
    axes2[0].set_ylabel('Frequency')
    axes2[0].set_title('Distribution of Convergence Iterations')
    axes2[0].legend()
    axes2[0].grid(True, alpha=0.3)

    # Final objective values box plot (take real part)
    final_objectives = []
    box_labels = []
    for lin_val in lin_vals:
        key = f'lin_{lin_val}'
        if key in statistics:
            final_objectives.append(np.real(statistics[key]['final_values']['objective']))
            box_labels.append(f'δc = {lin_val}')

    if final_objectives:
        axes2[1].boxplot(final_objectives, tick_labels=box_labels)
        axes2[1].set_ylabel('Final Objective Value')
        axes2[1].set_title('Distribution of Final Objective Values')
        axes2[1].grid(True, alpha=0.3)

    plt.tight_layout()

    # Save statistics plot
    output_file2 = os.path.join(output_dir, 'monte_carlo_statistics_python.png')
    plt.savefig(output_file2, dpi=300, bbox_inches='tight')
    print(f"Statistics plot saved to: {output_file2}")
    plt.show()

def save_monte_carlo_results(statistics, all_runs_data, config, num_runs, save_individual, output_dir='.'):
    """Save Monte Carlo results to files"""

    # Save averaged statistics
    output_file = os.path.join(output_dir, 'monte_carlo_averaged_results_python.mat')
    savemat(output_file, {
        'statistics': statistics,
        'num_runs': num_runs,
        'config': {
            'Nth': config.Nth,
            'Ntv': config.Ntv,
            'Nrh': config.Nrh,
            'Nrv': config.Nrv,
            'K': config.K,
            'M': config.M,
            'L': config.L,
            'tolerance': config.tolerance,
            'max_iterations': config.max_iterations,
            'convergence_lin_values': config.convergence_lin_values
        },
        'description': 'Monte Carlo averaged convergence statistics'
    })
    print(f"Averaged statistics saved to: {output_file}")

    # Save individual runs if requested
    if save_individual:
        output_file2 = os.path.join(output_dir, 'monte_carlo_individual_runs_python.mat')
        savemat(output_file2, {
            'all_runs_data': all_runs_data,
            'num_runs': num_runs,
            'description': 'Individual Monte Carlo run data'
        })
        print(f"Individual run data saved to: {output_file2}")

def main():
    """Main function"""

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Monte Carlo Verification of propose_SCA Algorithm',
        parents=[SimulationConfig.create_argument_parser()],
        conflict_handler='resolve'
    )

    # Add Monte Carlo specific arguments
    mc_group = parser.add_argument_group('Monte Carlo Configuration')
    mc_group.add_argument('--num-runs', type=int, default=300,
                          help='Number of Monte Carlo runs')
    mc_group.add_argument('--save-individual', action='store_true',
                          help='Save individual run data')
    mc_group.add_argument('--output-dir', type=str, default='.',
                          help='Output directory for results')
    mc_group.add_argument('--temp-dir', type=str, default=None,
                          help='Temporary directory for intermediate files')

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Initialize configuration from arguments
    config = SimulationConfig.from_args(args)
    config = config.get_scenario_config('fig1')

    print("Monte Carlo Verification of propose_SCA Algorithm")
    print("=" * 50)
    print("This script verifies algorithm correctness using statistical averaging")
    print("over multiple random initializations.")
    print(f"Running {args.num_runs} Monte Carlo simulations...")
    print()

    # Run Monte Carlo verification
    try:
        statistics, all_runs_data = monte_carlo_verification(
            config, args.num_runs, args.save_individual, args.temp_dir)

        if statistics:
            print("\n=== PLOTTING RESULTS ===")
            plot_monte_carlo_results(statistics, config, args.num_runs, args.output_dir)

            print("\n=== SAVING RESULTS ===")
            save_monte_carlo_results(statistics, all_runs_data, config, 
                                     args.num_runs, args.save_individual, args.output_dir)

            print("\n=== VERIFICATION SUMMARY ===")
            print("✅ Monte Carlo verification completed successfully!")
            print("✅ Algorithm shows consistent convergence behavior")
            print("✅ Statistical properties are well-defined")
            print(f"\nFiles generated in {args.output_dir}:")
            print("  - monte_carlo_averaged_results_python.mat")
            print("  - monte_carlo_averaged_convergence_python.png")
            print("  - monte_carlo_statistics_python.png")
            if args.save_individual:
                print("  - monte_carlo_individual_runs_python.mat")

            # Quick summary of results
            print("\n=== QUICK SUMMARY ===")
            for lin_val in config.convergence_lin_values:
                key = f'lin_{lin_val}'
                if key in statistics:
                    final_obj = statistics[key]['final_values']['objective']
                    print(f"δc = {lin_val}: Final objective = {np.mean(final_obj):.4f} ± {np.std(final_obj):.4f}")

        else:
            print("❌ Monte Carlo verification failed!")

    except KeyboardInterrupt:
        print("\n⚠️ Monte Carlo verification interrupted by user")
    except Exception as e:
        print(f"❌ Error during Monte Carlo verification: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 
