# statistics.py

"""Functions for statistical analysis of simulation results."""

import numpy as np
from collections import defaultdict
from optimization import bayesian_optimization, SSR_Score
from simulation import simulate_markov_model
from parameters import states, initial_state_distribution, observed_state_distribution_3days
from utils import remove_outliers, normalize_transition_probabilities

def run_simulation():
    """Run a single simulation and return the optimized probabilities and simulation history."""
    # Transition probabilities (with self-transition included)
    from parameters import Δ, VD, d0, VT, TD, Dr  # Import parameters

    transition_probabilities_hourly = {
        "S": np.array([1 - (VD + Δ), VD, 0, Δ, 0]),  # S to [S, E, I, R, D]
        "E": np.array([d0, 1 - (d0 + VT * TD + Δ), VT * TD + Δ, 0, 0]),
        "I": np.array([0, 0, 1 - (d0 + Dr), d0, Dr]),
        "R": np.array([0, 0, VD, 1 - (VD + Dr), Dr]),
        "D": np.array([0, 0, 0, 0, 1.0])
    }

    # Normalize the transition probabilities
    transition_probabilities_hourly = normalize_transition_probabilities(transition_probabilities_hourly)

    # Optimize the transition probabilities using Bayesian Optimization
    optimized_probabilities = bayesian_optimization(transition_probabilities_hourly, iterations=1000)

    # Simulate the model with optimized probabilities
    simulation_history = simulate_markov_model(optimized_probabilities, time_steps=72)

    return optimized_probabilities, simulation_history

def collect_and_compute_statistics(num_simulations=100):
    """Run multiple simulations and compute statistics.

    Args:
        num_simulations (int): Number of simulations to run.

    Returns:
        tuple: Mean probabilities, standard deviations, IQRs, mean SSR, std SSR, IQR SSR,
               mean final counts, all simulation histories.
    """
    probability_records = defaultdict(list)
    ssr_scores = []
    final_cell_counts = {state: [] for state in states}
    all_simulation_histories = []

    for _ in range(num_simulations):
        optimized_probabilities, simulation_history = run_simulation()

        # Collect probabilities
        for from_state in states:
            for to_state_index, to_state in enumerate(states):
                prob = optimized_probabilities[from_state][to_state_index]
                probability_records[(from_state, to_state)].append(prob)

        # Calculate SSR for this simulation
        final_distribution = simulation_history[-1]
        ssr = SSR_Score(final_distribution, observed_state_distribution_3days)
        ssr_scores.append(ssr)

        # Collect final cell counts
        for state in states:
            final_cell_counts[state].append(final_distribution[state])

        # Collect cell states over time
        all_simulation_histories.append(simulation_history)

    # Calculate statistics for probabilities
    mean_probabilities = {}
    std_probabilities = {}
    iqr_probabilities = {}

    for transition, probs in probability_records.items():
        # Remove outliers
        filtered_probs, IQR = remove_outliers(probs)
        mean_prob = np.mean(filtered_probs)
        std_prob = np.std(filtered_probs)
        mean_probabilities[transition] = mean_prob
        std_probabilities[transition] = std_prob
        iqr_probabilities[transition] = IQR

    # Calculate mean SSR
    mean_ssr = np.mean(ssr_scores)
    std_ssr = np.std(ssr_scores)
    iqr_ssr = np.percentile(ssr_scores, 75) - np.percentile(ssr_scores, 25)

    # Calculate mean final cell counts
    mean_final_counts = {state: np.mean(counts) for state, counts in final_cell_counts.items()}

    return (mean_probabilities, std_probabilities, iqr_probabilities,
            mean_ssr, std_ssr, iqr_ssr, mean_final_counts, all_simulation_histories)
