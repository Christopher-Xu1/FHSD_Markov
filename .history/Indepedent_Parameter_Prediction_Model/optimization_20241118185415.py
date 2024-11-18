"""Functions for optimizing transition probabilities.""

import numpy as np
from scipy.stats import dirichlet
from simulation import simulate_markov_model
from parameters import states, initial_state_distribution, observed_state_distribution_3days

def SSR_Score(predicted_distribution, observed_distribution):
    """Calculate the Sum of Squared Residuals (SSR) score.

    Args:
        predicted_distribution (dict): Predicted state distribution.
        observed_distribution (dict): Observed state distribution.

    Returns:
        float: SSR score.
    """
    score = sum((predicted_distribution[state] - observed_distribution[state]) ** 2 for state in predicted_distribution)
    return score

def bayesian_optimization(transition_probabilities, iterations=1000):
    """Optimize the transition probabilities using Bayesian Optimization.

    Args:
        transition_probabilities (dict): Initial transition probabilities.
        iterations (int): Number of optimization iterations.

    Returns:
        dict: Optimized transition probabilities.
    """
    best_probabilities = None
    best_score = float('inf')

    for _ in range(iterations):
        sampled_probabilities = {}
        for state, probs in transition_probabilities.items():
            if np.any(probs > 0):  # Only sample if there are non-zero probabilities
                non_zero_indices = probs > 0
                alpha = probs[non_zero_indices] + 1e-6  # Small value to avoid zero
                sampled_probs = dirichlet.rvs(alpha)[0]
                new_probs = np.zeros_like(probs)
                new_probs[non_zero_indices] = sampled_probs
                sampled_probabilities[state] = new_probs
            else:
                sampled_probabilities[state] = probs

        simulation_history = simulate_markov_model(sampled_probabilities, time_steps=72)  # 72 hours = 3 days
        final_distribution = simulation_history[-1]

        score = SSR_Score(final_distribution, observed_state_distribution_3days)
        if score < best_score:
            best_score = score
            best_probabilities = sampled_probabilities

    return best_probabilities
