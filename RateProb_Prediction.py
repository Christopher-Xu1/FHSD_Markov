import numpy as np
import matplotlib.pyplot as plt

# Initial and observed state distributions
initial_state_distribution = {"S": 5488, "E": 0, "I": 0, "R": 0, "D": 0}
observed_state_distribution_3days = {"S": 4956, "E": 14, "I": 13, "R": 150, "D": 355}

# Function to create transition probability matrix from rates
def create_transition_matrix(Δ, Dr, VD, d0, VT, TD):
    transition_probabilities_hourly = {
        "S": np.array([1 - (VD + Δ), VD, 0, Δ, 0]),
        "E": np.array([d0, 1 - (d0 + VT * TD + Δ), VT * TD + Δ, 0, 0]),
        "I": np.array([0, 0, 1 - (d0 + Dr), 0, d0]),
        "R": np.array([0, 0, VD, 1 - (VD + Dr), Dr]),
        "D": np.array([0, 0, 0, 0, 1.0])
    }
    return transition_probabilities_hourly

# Normalize non-zero probabilities
def normalize_probabilities(transition_probabilities):
    for key in transition_probabilities:
        total = transition_probabilities[key].sum()
        if total > 0:
            transition_probabilities[key] /= total  # Ensuring probabilities sum to 1
        else:
            transition_probabilities[key] = np.zeros_like(transition_probabilities[key])
    return transition_probabilities

# Function to check for invalid probabilities
def check_valid_probabilities(transition_probabilities):
    for state, probs in transition_probabilities.items():
        if not np.all(np.isfinite(probs)) or np.any(probs < 0) or np.any(probs > 1):
            return False
    return True

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

# Negative Log-Likelihood (NLL) function
def NLL_Score(predicted_distribution, observed_distribution):
    nll = 0
    for state in predicted_distribution:
        predicted_count = max(predicted_distribution[state], 1e-9)  # Avoid log(0)
        observed_count = observed_distribution[state]
        nll -= observed_count * np.log(predicted_count)
    return nll

# Bayesian Optimization with NLL tracking
def bayesian_optimization(transition_probabilities_func, initial_state_distribution, observed_distribution, iterations=10000):
    best_params = None
    best_score = float('inf')
    nll_scores = []

    for i in range(iterations):
        # Sample parameters randomly within some range
        Δ = np.random.uniform(0.001, 0.1)
        Dr = np.random.uniform(0.01, 0.1)
        VD = np.random.uniform(0.001, 0.01)
        d0 = np.random.uniform(0.1, 1.0)
        VT = np.random.uniform(1.0, 10.0)
        TD = np.random.uniform(0.01, 0.1)
        
        # Generate transition probabilities from sampled parameters
        transition_probabilities = transition_probabilities_func(Δ, Dr, VD, d0, VT, TD)
        transition_probabilities = normalize_probabilities(transition_probabilities)
        
        # Check if the transition probabilities are valid
        if not check_valid_probabilities(transition_probabilities):
            continue  # Skip this iteration if invalid probabilities are found

        # Simulate the model
        simulation_history = simulate_markov_model(transition_probabilities, initial_state_distribution, 72)  # 72 hours = 3 days
        final_distribution = simulation_history[-1]
    
        # Calculate NLL score
        score = NLL_Score(final_distribution, observed_distribution)
        nll_scores.append(score)

        if score < best_score:
            best_score = score
            best_params = (Δ, Dr, VD, d0, VT, TD)

        # Optional: print progress every 1000 iterations
        if (i + 1) % 1000 == 0:
            print(f"Iteration {i + 1}/{iterations}, Best NLL: {best_score}")
        

    # Plot the NLL scores over iterations
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(nll_scores)), nll_scores, label='NLL Score')
    plt.xlabel('Iteration')
    plt.ylabel('Negative Log-Likelihood (NLL)')
    plt.title('NLL Score Over Iterations')
    plt.grid(True)
    plt.show()

    return best_params

# Optimize the transition parameters using Bayesian Optimization
optimized_params = bayesian_optimization(create_transition_matrix, initial_state_distribution, observed_state_distribution_3days)

# Generate the optimized transition probabilities
optimized_transition_probabilities = create_transition_matrix(*optimized_params)
optimized_transition_probabilities = normalize_probabilities(optimized_transition_probabilities)

# Simulate the Markov model with optimized parameters
simulation_history = simulate_markov_model(optimized_transition_probabilities, initial_state_distribution, 72)

print("Optimized Parameters:")
print(f"Δ: {optimized_params[0]}")
print(f"Dr: {optimized_params[1]}")
print(f"VD: {optimized_params[2]}")
print(f"d0: {optimized_params[3]}")
print(f"VT: {optimized_params[4]}")
print(f"TD: {optimized_params[5]}")

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
