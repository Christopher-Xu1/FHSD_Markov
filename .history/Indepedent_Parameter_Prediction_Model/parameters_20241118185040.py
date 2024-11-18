# parameters.py
"""Module containing constants and parameters for the Markov model."""

# Biological Parameters
Δ = 0.01          # DUX4 syncytial diffusion rate
Dr = 1 / 20.1     # DUX4 target gene-induced death rate
VD = 0.00211      # Transcription rate
d0 = 0.246        # Degradation rate
VT = 6.41         # Translation rate
TD = 1 / 13       # mRNA half-life

# Initial and observed state distributions
initial_state_distribution = {"S": 5488, "E": 0, "I": 0, "R": 0, "D": 0}
observed_state_distribution_3days = {"S": 4956, "E": 14, "I": 13, "R": 150, "D": 355}

# States for transitions
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
