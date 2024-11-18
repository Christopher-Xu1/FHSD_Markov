"""Functions for simulating the Markov model."""

import numpy as np
from parameters import states, initial_state_distribution
from utils import normalize_transition_probabilities

def simulate_markov_model(transition_probabilities, time_steps):
    """Simulate the Markov model over time.

    Args:
        transition_probabilities (dict): Transition probabilities for each state.
        time_steps (int): Number of time steps to simulate.

    Returns:
        list: History of state distributions over time.
    """
    state_distribution = initial_state_distribution.copy()
    history = [state_distribution.copy()]

    for _ in range(time_steps):
        new_distribution = {state: 0 for state in states}
        for state, count in state_distribution.items():
            if count > 0:
                probs = transition_probabilities[state]
                transitions = np.random.multinomial(count, probs)
                for i, next_state in enumerate(states):
                    new_distribution[next_state] += transitions[i]
        state_distribution = new_distribution.copy()
        history.append(state_distribution.copy())

    return history
