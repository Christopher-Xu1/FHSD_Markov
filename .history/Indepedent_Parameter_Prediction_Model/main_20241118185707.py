from statistics import collect_and_compute_statistics
from display import print_transition_table
from plotting import plot_mean_cell_states, plot_final_distribution
from parameters import states

def main():
    """Main function to run simulations and display results."""
    num_simulations = 100  # Number of simulations to run

    # Collect and compute statistics over simulations
    (mean_probabilities, std_probabilities, iqr_probabilities,
     mean_ssr, std_ssr, iqr_ssr,
     mean_final_counts, all_simulation_histories) = collect_and_compute_statistics(num_simulations)

    # Print the transition table
    print_transition_table(mean_probabilities, std_probabilities, iqr_probabilities)

    # Print the mean SSR
    print(f"\nMean SSR over simulations: {mean_ssr:.5f}")
    print(f"Standard Deviation of SSR: {std_ssr:.5f}")
    print(f"IQR of SSR: {iqr_ssr:.5f}")

    # Print the mean final cell counts
    print("\nMean Final Cell Counts over simulations:")
    for state in states:
        print(f"{state}: {mean_final_counts[state]:.2f}")

    # Plot the mean cell states over time
    plot_mean_cell_states(all_simulation_histories)

    # Plot the final distributions for comparison
    from parameters import observed_state_distribution_3days
    plot_final_distribution(mean_final_counts, observed_state_distribution_3days)

if __name__ == "__main__":
    main()
