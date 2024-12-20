import numpy as np
from scipy.stats import dirichlet, chisquare
import matplotlib.pyplot as plt
import random
# Data retrieval
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix
from tabulate import tabulate


# # Define the sample file mappings (local paths)
# samples = {
#     "FSHD1.1": "/Users/chris/iGEM_modeling/scRNAseqData/GSM3487556_FSHD1.1.txt",
#     "FSHD1.2": "/Users/chris/iGEM_modeling/scRNAseqData/GSM3487557_FSHD1.2.txt"
#     }

# # Initialize an empty dictionary to store the AnnData objects
# adatas = {}

# # Read and process data for each sample
# for sample_id, filepath in samples.items():
#     try:
#         # Read the text file into a DataFrame 
#         sample_data = pd.read_csv(filepath, sep="\t", index_col=0)
        
#         # Convert the DataFrame to a sparse matrix
#         sample_matrix = csr_matrix(sample_data.values)
#         # Create an AnnData object
#         sample_adata = ad.AnnData(sample_matrix)

#         # Set the observation (cell) names and variable (gene) names
#         sample_adata.var_names_make_unique()  # Ensure gene names are unique
#         sample_adata.obs_names = sample_data.index.tolist()
#         sample_adata.var_names = sample_data.columns.tolist()
        
#         # Store the AnnData object in the dictionary
#         adatas[sample_id] = sample_adata
        
#         print(f"Successfully read data for {sample_id}")
#     except Exception as e:
#         print(f"Failed to read data for {sample_id}: {e}")

# # The adatas dictionary now contains AnnData objects for each sample
# for sample_id, adata in adatas.items():
#     print(f"Data for {sample_id}:")
#     print(adata)  # Print a summary of each AnnData object
#     print(adata.X)  # Print the data matrix

# # Optional: Print the observation and variable names for the first sample
# if adatas:
#     sample_id = list(adatas.keys())[0]
#     adata = adatas[sample_id]
#     print("Observation (cell) names:", adata.obs_names[:10])
#     print("Variable (gene) names:", adata.var_names[:10])

# adata.layers["log_transformed"] = np.log1p(adata.X)
# adata

# adata.to_df(layer="log_transformed")














# Initial and observed state distributions
initial_state_distribution = {"S": 5488, "E": 0, "I": 0, "R": 0, "D": 0}
observed_state_distribution_3days = {"S": 4956, "E": 14, "I": 13, "R": 150, "D": 355}
import numpy as np
from scipy.stats import dirichlet
import matplotlib.pyplot as plt

# Initial and observed state distributions
initial_state_distribution = {"S": 5488, "E": 0, "I": 0, "R": 0, "D": 0}
observed_state_distribution_3days = {"S": 4956, "E": 14, "I": 13, "R": 150, "D": 355}

import random
import numpy as np

# Parameters
Δ = 0.01  # DUX4 syncytial diffusion rate
Dr = 1/20.1  # DUX4 target gene-induced death rate
VD = 0.00211  # Transcription rate
d0 = 0.246  # Degradation rate
VT = 6.41  # Translation rate
TD = 1/13  # mRNA half-life

# Transition probabilities (with self-transition included)
transition_probabilities_hourly = {
    "S": np.array([1 - (VD + Δ), VD, 0, Δ, 0]),  # S to [S, E, I, R, D]
    "E": np.array([d0, 1 - (d0 + VT * TD + Δ), VT * TD + Δ, 0, 0]),
    "I": np.array([0, 0, 1 - (d0 + Dr), d0, Dr]),
    "R": np.array([0, 0, VD, 1 - (VD + Dr), Dr]),
    "D": np.array([0, 0, 0, 0, 1.0])
}

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

# Normalize non-zero probabilities
for key in transition_probabilities_hourly:
    transition_probabilities_hourly[key] /= transition_probabilities_hourly[key].sum()  # Ensuring probabilities sum to 1
# Function to simulate Markov model over time
def simulate_markov_model(transition_probabilities, initial_state_distribution, time_steps):
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

# Bayesian Optimization
def bayesian_optimization(transition_probabilities, initial_state_distribution, observed_distribution, iterations=10000):
    best_probabilities = None
    best_score = float('inf')

    for _ in range(iterations):
        sampled_probabilities = {}
        for state, probs in transition_probabilities.items():
            if np.any(probs > 0):  # Only sample from Dirichlet if there are non-zero probabilities
                non_zero_indices = probs > 0
                sampled_probs = dirichlet.rvs(probs[non_zero_indices])[0]
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

def SSR_Score(predicted_distribution, observed_distribution):
    score = sum((predicted_distribution[state] - observed_distribution[state]) ** 2 for state in predicted_distribution)
    return score

# Optimize the transition probabilities using Bayesian Optimization
optimized_probabilities = bayesian_optimization(transition_probabilities_hourly, initial_state_distribution, observed_state_distribution_3days)
simulation_history = simulate_markov_model(optimized_probabilities, initial_state_distribution, 72)


rounded_probabilities = {
    state: [round(prob, 5) for prob in probabilities]
    for state, probabilities in optimized_probabilities.items()
}

# Convert to a table format for neat display
table = []
for state, probabilities in rounded_probabilities.items():
    table.append([state] + probabilities)

# Print the table with proper alignment
headers = ['State'] + ['S']+ ['E']+ ['I'] + ['R'] + ['D']
print(tabulate(table, headers=headers, tablefmt="grid"))



# Prepare the table data with formatted cells
table_data = []
header = ["From\\To"] + states

for from_state in states:
    row = [from_state]
    for to_state in states:
        prob = optimized_probabilities[from_state][states.index(to_state)]
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

# Print the optimized transition probabilities table
print("Optimized Transition Probabilities:")
print_table(header, table_data)




# Print the history of state distributions at each hour
# for t, distribution in enumerate(simulation_history):
#     print(f"Hour {t} Cell States: {distribution}")


# Validate the final state distribution with observed data at 72 hours (3 days)
final_distribution = simulation_history[-1]
observed_distribution = observed_state_distribution_3days

print("\nFinal simulated distribution vs. observed distribution (72 hours):")
for state in final_distribution:
    print(f"{state}: Simulated={final_distribution[state]}, Observed={observed_distribution[state]}")





# Plot the cell states over time
time_points = range(len(simulation_history))
states = list(simulation_history[0].keys())
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

# Plot the final distributions for comparison
states = list(final_distribution.keys())
simulated_counts = [final_distribution[state] for state in states]
observed_counts = [observed_distribution[state] for state in states]

x = np.arange(len(states))
width = 0.35

fig, ax = plt.subplots()
bars1 = ax.bar(x - width/2, simulated_counts, width, label='Simulated')
bars2 = ax.bar(x + width/2, observed_counts, width, label='Observed')

ax.set_ylabel('Counts')
ax.set_title('Final State Distribution at 72 Hours')
ax.set_xticks(x)
ax.set_xticklabels(states)
ax.legend()

plt.show()
