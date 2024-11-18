"""Functions for plotting simulation results."""

import matplotlib.pyplot as plt
from parameters import states

def plot_mean_cell_states(all_simulation_histories):
    """Plot the mean cell states over time across all simulations.

    Args:
        all_simulation_histories (list): List of simulation histories.
    """
    num_simulations = len(all_simulation_histories)
    time_steps = len(all_simulation_histories[0])
    mean_state_counts = {state: [0] * time_steps for state in states}

    for simulation_history in all_simulation_histories:
        for t, distribution in enumerate(simulation_history):
            for state in states:
                mean_state_counts[state][t] += distribution[state]

    # Calculate the mean over simulations
    for state in states:
        mean_state_counts[state] = [count / num_simulations for count in mean_state_counts[state]]

    # Plotting
    time_points = range(time_steps)
    plt.figure(figsize=(10, 6))

    for state in states:
        plt.plot(time_points, mean_state_counts[state], label=state)

    plt.xlabel('Time (hours)')
    plt.ylabel('Average Number of Cells')
    plt.title('Mean Cell States Over Time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_final_distribution(mean_final_counts, observed_distribution):
    """Plot the mean simulated and observed distributions.

    Args:
        mean_final_counts (dict): Mean final counts from simulations.
        observed_distribution (dict): Observed state distribution.
    """
    states_list = list(mean_final_counts.keys())
    simulated_counts = [mean_final_counts[state] for state in states_list]
    observed_counts = [observed_distribution[state] for state in states_list]

    x = range(len(states_list))
    width = 0.35

    fig, ax = plt.subplots()
    ax.bar([xi - width / 2 for xi in x], simulated_counts, width, label='Mean Simulated')
    ax.bar([xi + width / 2 for xi in x], observed_counts, width, label='Observed')

    ax.set_ylabel('Counts')
    ax.set_title('Mean Final State Distribution at 72 Hours')
    ax.set_xticks(x)
    ax.set_xticklabels(states_list)
    ax.legend()

    plt.show()
