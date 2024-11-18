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


Δ = 0.01  # DUX4 syncytial diffusion rate
Dr = 1/20.1  # DUX4 target gene-induced death rate
VD = 0.00211  # Transcription rate
d0 = 0.246  # Degradation rate
VT = 6.41 # Translation rate
TD = 1/13  # mRNA half-life
xs = random.uniform(0.1, 1)
xe = random.uniform(0.1, 1)
xi = random.uniform(0.1, 1)
xr = random.uniform(0.1, 1)

# Hourly transition probabilities (with self-transition included)
transition_probabilities_hourly = {
    "S": np.array([1 - (VD + Δ), VD, 0, Δ, 0]),  # Δ influences transition from S to R
    "E": np.array([d0, 1 - (d0 + VT * TD + Δ), VT * TD + Δ, 0, 0]),
    "I": np.array([0, 0, 1 - (d0 + Dr), d0, Dr]),
    "R": np.array([0, 0, VD, 1 - (VD + Dr), Dr]),
    "D": np.array([0, 0, 0, 0, 1.0])

}

import random
import numpy as np

# Parameters
Δ = 0.01  # DUX4 syncytial diffusion rate
Dr = 1/20.1  # DUX4 target gene-induced death rate
VD = 0.00211  # Transcription rate
d0 = 0.246  # Degradation rate
VT = 6.41  # Translation rate
TD = 1/13  # mRNA half-life

# Hourly transition probabilities (with self-transition included)
transition_probabilities_hourly = {
    "S": np.array([1 - (VD + Δ), VD, 0, Δ, 0]),  # Δ influences transition from S to R
    "E": np.array([d0, 1 - (d0 + VT * TD + Δ), VT * TD + Δ, 0, 0]),
    "I": np.array([0, 0, 1 - (d0 + Dr), d0, Dr]),
    "R": np.array([0, 0, VD, 1 - (VD + Dr), Dr]),
    "D": np.array([0, 0, 0, 0, 1.0])
}

# Function to print transition probabilities with corresponding rates
def print_transition_table(transition_probabilities):
    rate_labels = {
        "S": [f"Δ = {Δ}", f"VD = {VD}", "None", f"Δ = {Δ}", "None"],
        "E": [f"d0 = {d0}", f"Δ = {Δ}", f"VT * TD + Δ = {VT * TD + Δ}", "None", "None"],
        "I": ["0", "0", f"Dr = {Dr}", f"d0 = {d0}", f"Dr = {Dr}"],
        "R": ["0", "None", f"VD = {VD}", f"Dr = {Dr}", f"Dr = {Dr}"],
        "D": ["0", "0", "0", "0", "1.0"]
    }

    # Print the header
    print(f"{'State':<10}", end="")  # Left-aligned header for State
    for i in range(len(next(iter(transition_probabilities.values())))):
        print(f"Prob {i+1:<5}", end="")  # Left-aligned header for probabilities
    print()  # New line after the header

    # Print each state, its probabilities, and corresponding rates
    for state, probabilities in transition_probabilities.items():
        print(f"{state:<10}", end="")  # Left-aligned state
        for prob, rate in zip(probabilities, rate_labels[state]):
            print(f"{prob:<6.5f} ({rate})", end="")  # Left-aligned rounded probability with rate
        print()  # New line after each state

# Print the transition table with labels
print_transition_table(transition_probabilities_hourly)


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



# Print the history of state distributions at each hour
# for t, distribution in enumerate(simulation_history):
#     print(f"Hour {t} Cell States: {distribution}")

# Plot the cell states over time
time_points = range(len(simulation_history))
states = list(simulation_history[0].keys())

plt.figure(figsize=(10, 6))

#Print States Over Time
for state in states:
    state_counts = [distribution[state] for distribution in simulation_history]
    plt.plot(time_points, state_counts, label=state)

plt.xlabel('Time (hours)')
plt.ylabel('Number of Cells')
plt.title('Cell States Over Time')
plt.legend()
plt.grid(True)
plt.show()

# Validate the final state distribution with observed data at 72 hours (3 days)
final_distribution = simulation_history[-1]
observed_distribution = observed_state_distribution_3days

print("\nFinal simulated distribution vs. observed distribution (72 hours):")
for state in final_distribution:
    print(f"{state}: Simulated={final_distribution[state]}, Observed={observed_distribution[state]}")

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
