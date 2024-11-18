# utils.py

"""Utility functions for the Markov model simulation."""

import numpy as np

def normalize_transition_probabilities(transition_probabilities):
    """Normalize non-zero probabilities to ensure they sum to 1."""
    for key in transition_probabilities:
        total = transition_probabilities[key].sum()
        if total > 0:
            transition_probabilities[key] /= total
    return transition_probabilities

def remove_outliers(data):
    """Remove outliers from a list of data using the IQR method."""
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    filtered_data = [x for x in data if lower_bound <= x <= upper_bound]
    return filtered_data, IQR
