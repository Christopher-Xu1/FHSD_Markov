import numpy as np
from scipy.stats import dirichlet
import matplotlib.pyplot as plt
import random
from collections import defaultdict
import statistics

# Parameters (constants)
Δ = 0.01  # DUX4 syncytial diffusion rate
Dr = 1 / 20.1  # DUX4 target gene-induced death rate
VD = 0.00211  # Transcription rate
d0 = 0.246  # Degradation rate
VT = 6.41  # Translation rate
TD = 1 / 13  # mRNA half-life

# Initial and observed state distributions
initial_state_distribution = {"S": 5488, "E": 0, "I": 0, "R": 0, "D": 0}
observed_state_distribution_3days = {"S": 4956, "E": 14, "I": 13, "R": 150, "D": 355}

# Define the states for transitions
states = ["S", "E", "I", "R", "D"]

# Mapping of transitions to their corresponding concise rate labels
transition_rates = {
    ("S", "S"): "",  # Self-transition
    ("S", "E"): "VD",
    ("S", "I"): "",
    ("S", "R"): "Δ",
    ("S", "D"): "",

    ("E", "S"): "d0",
    ("E", "E"): "",
    ("E", "I"): "VT*TD+Δ",
    ("E", "R"): "",
    ("E", "D"): "",

    ("I", "S"): "",
    ("I", "E"): "",
    ("I", "I"): "",
    ("I", "R"): "d0",
    ("I", "D"): "Dr",

    ("R", "S"): "",
    ("R", "E"): "",
    ("R", "I"): "VD",
    ("R", "R"): "",
    ("R", "D"): "Dr",

    ("D", "S"): "",
    ("D", "E"): "",
    ("D", "I"): "",
    ("D", "R"): "",
    ("D", "D"): ""
}

def normalize_transition_probabilities(transition_probabilities):
    """Normalize non-zero probabilities to ensure they sum to 1."""
    for key in transition_probabilities:
        total = transition_probabilities[key].sum()
        if total > 0:
            transition_probabilities[key] /= total
    return transition_probabilities

def simulate_markov_model(transition_probabilities, initial_state_distribution, time_steps):
    """Simulate the Markov model over time."""
    state_distribution = initial_state_distribution.copy()
    states = list(initial_state_distribution.keys())
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

def SSR_Score(predicted_distribution, observed_distribution):
    """Calculate the Sum of Squared Residuals (SSR) score."""
    score = sum((predicted_distribution[state] - observed_distribution[state]) ** 2 for state in predicted_distribution)
    return score

def bayesian_optimization(transition_probabilities, initial_state_distribution, observed_distribution, iterations=10000):
    """Optimize the transition probabilities using Bayesian Optimization."""
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

        simulation_history = simulate_markov_model(sampled_probabilities, initial_state_distribution, 72)  # 72 hours = 3 days
        final_distribution = simulation_history[-1]

        score = SSR_Score(final_distribution, observed_distribution)
        if score < best_score:
            best_score = score
            best_probabilities = sampled_probabilities

    return best_probabilities

def run_simulation():
    """Run a single simulation and return the optimized probabilities and simulation history."""
    # Transition probabilities (with self-transition included)
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
    optimized_probabilities = bayesian_optimization(transition_probabilities_hourly, initial_state_distribution, observed_state_distribution_3days, iterations=1000)

    # Simulate the model with optimized probabilities
    simulation_history = simulate_markov_model(optimized_probabilities, initial_state_distribution, 72)

    return optimized_probabilities, simulation_history

def collect_and_compute_median_probabilities(num_simulations=100):
    """Run multiple simulations and compute the median of the optimized probabilities."""
    # Initialize a nested dictionary to store probabilities
    probability_records = defaultdict(list)

    for _ in range(num_simulations):
        optimized_probabilities, _ = run_simulation()
        for from_state in states:
            for to_state_index, to_state in enumerate(states):
                prob = optimized_probabilities[from_state][to_state_index]
                probability_records[(from_state, to_state)].append(prob)

    # Calculate median probabilities
    median_probabilities = {}
    for transition, probs in probability_records.items():
        median_prob = statistics.median(probs)
        median_probabilities[transition] = median_prob

    # Prepare median transition probabilities
    median_transition_probabilities = {}
    for from_state in states:
        median_transition_probabilities[from_state] = np.array([
            median_probabilities[(from_state, to_state)] for to_state in states
        ])

    return median_transition_probabilities

def print_transition_table(transition_probabilities):
    """Print the transition probabilities in a formatted table."""
    # Round the probabilities for display
    rounded_probabilities = {
        state: [round(prob, 5) for prob in probabilities]
        for state, probabilities in transition_probabilities.items()
    }

    # Prepare the table data with formatted cells
    table_data = []
    header = ["From\\To"] + states

    for from_state in states:
        row = [from_state]
        for to_state in states:
            prob = rounded_probabilities[from_state][states.index(to_state)]
            prob_str = f"{prob:.5f}"
            rate_label = transition_rates.get((from_state, to_state), "")
            if rate_label:
                cell = f"{prob_str}\n({rate_label})"
            else:
                cell = prob_str
            row.append(cell)
        table_data.append(row)

    # Function to calculate column widths
    def calculate_column_widths(table, headers):
        column_widths = [len(h) for h in headers]
        for row in table:
            for i, cell in enumerate(row):
                cell_lines = cell.split('\n')
                max_line_length = max(len(line) for line in cell_lines)
                if max_line_length > column_widths[i]:
                    column_widths[i] = max_line_length
        return column_widths

    # Function to print the table with alignment
    def print_table(headers, table):
        column_widths = calculate_column_widths(table, headers)

        def print_separator():
            line = '+'
            for width in column_widths:
                line += '-' * (width + 2) + '+'
            print(line)

        def print_row(row):
            max_lines = max(cell.count('\n') + 1 for cell in row)
            cell_lines_list = [cell.split('\n') for cell in row]
            for i in range(max_lines):
                line = '|'
                for j, cell_lines in enumerate(cell_lines_list):
                    cell_line = cell_lines[i] if i < len(cell_lines) else ''
                    line += ' ' + cell_line.center(column_widths[j]) + ' |'
                print(line)

        print_separator()
        print_row(headers)
        print_separator()
        for row in table:
            print_row(row)
            print_separator()

    # Print the transition probabilities table
    print("Median Transition Probabilities over Simulations:")
    print_table(header, table_data)

def main():
    """Main function to run simulations and display results."""
    num_simulations = 100  # Number of simulations to run

    # Collect and compute median probabilities over simulations
    median_transition_probabilities = collect_and_compute_median_probabilities(num_simulations)

    # Simulate the model using median probabilities
    simulation_history = simulate_markov_model(median_transition_probabilities, initial_state_distribution, 72)

    # Print the transition table
    print_transition_table(median_transition_probabilities)

    # Print the final simulated distribution vs. observed distribution
    final_distribution = simulation_history[-1]
    observed_distribution = observed_state_distribution_3days

    print("\nFinal simulated distribution vs. observed distribution (72 hours):")
    for state in final_distribution:
        print(f"{state}: Simulated={final_distribution[state]}, Observed={observed_distribution[state]}")

    # Plot the cell states over time
    plot_simulation_history(simulation_history)

    # Plot the final distributions for comparison
    plot_final_distribution(final_distribution, observed_distribution)

def plot_simulation_history(simulation_history):
    """Plot the cell states over time."""
    time_points = range(len(simulation_history))
    plt.figure(figsize=(10, 6))

    for state in states:
        state_counts = [distribution[state] for distribution in simulation_history]
        plt.plot(time_points, state_counts, label=state)

    plt.xlabel('Time (hours)')
    plt.ylabel('Number of Cells')
    plt.title('Cell States Over Time')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_final_distribution(final_distribution, observed_distribution):
    """Plot the final simulated and observed distributions."""
    states_list = list(final_distribution.keys())
    simulated_counts = [final_distribution[state] for state in states_list]
    observed_counts = [observed_distribution[state] for state in states_list]

    x = np.arange(len(states_list))
    width = 0.35

    fig, ax = plt.subplots()
    bars1 = ax.bar(x - width / 2, simulated_counts, width, label='Simulated')
    bars2 = ax.bar(x + width / 2, observed_counts, width, label='Observed')

    ax.set_ylabel('Counts')
    ax.set_title('Final State Distribution at 72 Hours')
    ax.set_xticks(x)
    ax.set_xticklabels(states_list)
    ax.legend()

    plt.show()

if __name__ == "__main__":
    main()
