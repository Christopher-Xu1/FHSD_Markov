"""Functions for displaying outputs, including printing tables."""

from parameters import states, transition_rates

def print_transition_table(mean_probabilities, std_probabilities, iqr_probabilities):
    """Print the mean transition probabilities in a formatted table.

    Args:
        mean_probabilities (dict): Mean transition probabilities.
        std_probabilities (dict): Standard deviations of the probabilities.
        iqr_probabilities (dict): Interquartile ranges of the probabilities.
    """
    # Prepare mean transition probabilities
    mean_transition_probabilities = {}
    for from_state in states:
        mean_transition_probabilities[from_state] = [
            mean_probabilities[(from_state, to_state)] for to_state in states
        ]

    # Prepare the table data with formatted cells
    table_data = []
    header = ["From\\To"] + states

    for from_state in states:
        row = [from_state]
        for to_state_index, to_state in enumerate(states):
            prob_mean = mean_transition_probabilities[from_state][to_state_index]
            prob_std = std_probabilities[(from_state, to_state)]
            prob_iqr = iqr_probabilities[(from_state, to_state)]
            prob_str = f"{prob_mean:.5f}\n±{prob_std:.5f}\nIQR={prob_iqr:.5f}"
            rate_label = transition_rates.get((from_state, to_state), "")
            if rate_label:
                cell = f"{prob_str}\n({rate_label})"
            else:
                cell = prob_str
            row.append(cell)
        table_data.append(row)

    # Print the table
    print_table(header, table_data)

def print_table(headers, table):
    """Prints a formatted table.

    Args:
        headers (list): List of header strings.
        table (list): List of rows, where each row is a list of cell strings.
    """
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

    print("Mean Transition Probabilities over Simulations (Outliers Removed):")
    print_separator()
    print_row(headers)
    print_separator()
    for row in table:
        print_row(row)
        print_separator()

def calculate_column_widths(table, headers):
    """Calculates the width of each column based on the table data.

    Args:
        table (list): Table data.
        headers (list): Header strings.

    Returns:
        list: Column widths.
    """
    column_widths = [len(h) for h in headers]
    for row in table:
        for i, cell in enumerate(row):
            cell_lines = cell.split('\n')
            max_line_length = max(len(line) for line in cell_lines)
            if max_line_length > column_widths[i]:
                column_widths[i] = max_line_length
    return column_widths
