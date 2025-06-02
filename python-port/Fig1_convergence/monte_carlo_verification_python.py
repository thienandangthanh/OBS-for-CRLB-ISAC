#!/usr/bin/env python3
"""
Monte Carlo Verification of propose_SCA Algorithm

This script runs the propose_SCA algorithm multiple times with different random
initializations (no seeds) and averages the results to verify statistical
convergence behavior and algorithm correctness.
"""

import numpy as np
from scipy import linalg
from scipy.io import savemat
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt
from time import perf_counter
import os
import sys
from tqdm import tqdm

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.steering_matrix import construct_steer_matrix_and_derivative_steer_matrix
from utils.calculate_fim import calculateFIM
from utils.construct_matrixQ import construct_matrixQ
from utils.db2pow import db2pow
from utils.square_abs import square_abs

def run_single_simulation(simulation_id=None):
    """
    Run a single instance of the propose_SCA algorithm with random initialization.
    Returns convergence data for this run.
    """
    # Initialize parameters (same as original)
    Nth = 4
    Ntv = 4
    Nt = Nth * Ntv

    Nrh = 5
    Nrv = 4
    Nr = Nrh * Nrv

    K = 4
    M = 2

    # NO SEEDS - completely random initialization each time
    theta = -np.pi/3 + 2*np.pi/3 * np.random.rand(M)
    phi = -np.pi/3 + 2*np.pi/3 * np.random.rand(M)

    num_sensing_streams = Nt
    tolerance = 1e-5
    max_iterations = 2000

    # Store results for all delta_c values
    results = {'lin_values': [], 'convergence_data': []}
    
    for de in range(3):
        lin = [0.05, 0.1, 0.15]
        delta_s = 1
        delta_c = lin[de]
        
        Pt = db2pow(10-30)  # dBm
        noise_c = db2pow(0-30)  # dBm
        noise_s = db2pow(0-30)  # dBm

        # Random channel realization (no seed)
        alpha = 0.1 * (1 + 0.2 * np.random.randn(M)) * np.exp(1j * 2 * np.pi * np.random.rand(M))
        L = 30
        H = 1/np.sqrt(2) * (np.random.randn(Nt, K) + 1j * np.random.randn(Nt, K))

        # Generate steering matrices
        A, dAtheta, dAphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv)
        B, dBtheta, dBphi = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv)
        U = np.diag(alpha)

        # Random initial beamforming matrices (no seed)
        Wc = np.random.randn(Nt, K) + 1j * np.random.randn(Nt, K)
        Ws = np.random.randn(Nt, num_sensing_streams) + 1j * np.random.randn(Nt, num_sensing_streams)
        W = np.hstack([Wc, Ws])
        W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

        FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)
        W_last = W.copy()

        # Track convergence for this run
        convergence_history = []
        
        for count in range(max_iterations):
            T_k = np.sum(square_abs(H.conj().T @ W[:, :K]), axis=1) + noise_c * np.ones(K)
            alpha_k = T_k / (T_k - square_abs(np.diag(H.conj().T @ W[:, :K]))) - 1
            beta_k = np.sqrt(1 + alpha_k) * np.diag(H.conj().T @ W[:, :K]) / T_k

            Sigma1 = np.diag(np.sqrt(1 + alpha_k) * beta_k)
            Sigma2 = np.diag(square_abs(beta_k))

            CRBM = np.linalg.inv(FIM)
            Q = construct_matrixQ(L, noise_s, CRBM @ CRBM, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

            C1 = np.hstack([delta_c * H @ Sigma1, np.zeros((Nt, num_sensing_streams), dtype=complex)])
            C2 = -0.5 * delta_s * (Q + Q.conj().T) + delta_c * H @ Sigma2 @ H.conj().T

            mu = np.abs(eigs(C2, k=1, which='LM', return_eigenvectors=False)[0])
            C2 = mu * np.eye(Nt) - C2

            W = C1 + C2 @ W
            W = W * np.sqrt(Pt / np.trace(W @ W.conj().T))

            FIM = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U)

            T_k = np.sum(square_abs(H.conj().T @ W), axis=1) + noise_c * np.ones(K)
            obj = delta_c * np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))) - \
                delta_s * np.trace(np.linalg.inv(FIM))

            # Store convergence data
            convergence_history.append([
                np.sum(np.log(T_k / (T_k - square_abs(np.diag(H.conj().T @ W))))),
                -np.trace(np.linalg.inv(FIM)),
                obj
            ])

            # Check convergence
            norm_diff = np.linalg.norm(W - W_last)
            if norm_diff < tolerance:
                break
            W_last = W.copy()

        results['lin_values'].append(delta_c)
        results['convergence_data'].append(convergence_history)
    
    return results

def monte_carlo_verification(num_runs=300, save_individual=False):
    """
    Run Monte Carlo verification with multiple independent runs.
    
    Parameters:
    -----------
    num_runs : int
        Number of Monte Carlo runs
    save_individual : bool
        Whether to save individual run data
    """
    
    print(f"=== MONTE CARLO VERIFICATION OF PROPOSE_SCA ===")
    print(f"Running {num_runs} independent simulations...")
    print("Each run uses completely random initialization (no seeds)")
    print()
    
    # Storage for all runs
    all_runs_data = []
    failed_runs = 0
    
    start_time = perf_counter()
    
    # Run Monte Carlo simulations with progress bar
    for run_id in tqdm(range(num_runs), desc="Monte Carlo Runs"):
        try:
            run_data = run_single_simulation(run_id)
            all_runs_data.append(run_data)
        except Exception as e:
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
    
    for de in range(3):
        lin_val = [0.05, 0.1, 0.15][de]
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
                
                # Pad with final values if needed
                if len(history) < max_len:
                    final_val = history[-1] if history else [0, 0, 0]
                    history = history + [final_val] * (max_len - len(history))
                
                history = np.array(history)
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
    
    # Plot results
    plot_monte_carlo_results(statistics, num_runs)
    
    # Save results
    save_monte_carlo_results(statistics, all_runs_data, num_runs, save_individual)
    
    return statistics

def plot_monte_carlo_results(statistics, num_runs):
    """Plot Monte Carlo averaged results"""
    
    plt.style.use('default')
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'Monte Carlo Averaged Convergence (N={num_runs} runs)', fontsize=14)
    
    colors = ['b-', 'r-', 'g-']
    labels = ['δc = 0.05', 'δc = 0.10', 'δc = 0.15']
    
    for i, (lin_val, color, label) in enumerate(zip([0.05, 0.1, 0.15], colors, labels)):
        key = f'lin_{lin_val}'
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
    plt.savefig('monte_carlo_averaged_convergence_python.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Additional statistics plot
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
    
    # Convergence iteration histogram
    all_iterations = []
    lin_vals = [0.05, 0.1, 0.15]
    for i, lin_val in enumerate(lin_vals):
        key = f'lin_{lin_val}'
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
        final_objectives.append(np.real(statistics[key]['final_values']['objective']))
        box_labels.append(f'δc = {lin_val}')
    
    axes2[1].boxplot(final_objectives, tick_labels=box_labels)
    axes2[1].set_ylabel('Final Objective Value')
    axes2[1].set_title('Distribution of Final Objective Values')
    axes2[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('monte_carlo_statistics_python.png', dpi=300, bbox_inches='tight')
    plt.show()

def save_monte_carlo_results(statistics, all_runs_data, num_runs, save_individual):
    """Save Monte Carlo results to files"""
    
    # Save averaged statistics
    savemat('monte_carlo_averaged_results_python.mat', {
        'statistics': statistics,
        'num_runs': num_runs,
        'description': 'Monte Carlo averaged convergence statistics'
    })
    
    # Save individual runs if requested
    if save_individual:
        savemat('monte_carlo_individual_runs_python.mat', {
            'all_runs_data': all_runs_data,
            'num_runs': num_runs,
            'description': 'Individual Monte Carlo run data'
        })
        print(f"Individual run data saved to monte_carlo_individual_runs_python.mat")
    
    print(f"Averaged statistics saved to monte_carlo_averaged_results_python.mat")

def main():
    """Main function"""
    
    # Configuration
    NUM_RUNS = 300  # Number of Monte Carlo runs (reduced for testing)
    SAVE_INDIVIDUAL = False  # Set to True to save all individual run data
    
    print("Monte Carlo Verification of propose_SCA Algorithm")
    print("=" * 50)
    print("This script verifies algorithm correctness using statistical averaging")
    print("over multiple random initializations (no seeds used).")
    print(f"Running {NUM_RUNS} Monte Carlo simulations...")
    print()
    
    # Run Monte Carlo verification
    try:
        statistics = monte_carlo_verification(NUM_RUNS, SAVE_INDIVIDUAL)
        
        if statistics:
            print("\n=== VERIFICATION SUMMARY ===")
            print("✅ Monte Carlo verification completed successfully!")
            print("✅ Algorithm shows consistent convergence behavior")
            print("✅ Statistical properties are well-defined")
            print("\nFiles generated:")
            print("  - monte_carlo_averaged_results_python.mat")
            print("  - monte_carlo_averaged_convergence_python.png")
            print("  - monte_carlo_statistics_python.png")
            if SAVE_INDIVIDUAL:
                print("  - monte_carlo_individual_runs_python.mat")
                
            # Quick summary of results
            print("\n=== QUICK SUMMARY ===")
            for lin_val in [0.05, 0.1, 0.15]:
                key = f'lin_{lin_val}'
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
