# CREATE THE DEACTIVATE HAIRPINS

import os
import sys

if __name__=="__main__":
    user_file = sys.argv[1]
    user_file2 = sys.argv[2]

input_file = user_file #"temp_classes_table.tsv"
output_file = user_file + "_temp_activated_hairpins.txt"

# Open the input file for reading and output file for writing
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    # Iterate through each line in the file
    for line in infile:
        # Check if the line contains the desired status
        if "\tactivated" in line:
            # Split the line into columns by tab
            columns = line.strip().split("\t")
            # Extract the first and last columns
            first_column = columns[0]
            last_column = columns[-1]
            # Write the extracted columns to the output file
            outfile.write(f"{first_column}\t{last_column}\n")


# DO THE FISHER TEST

from collections import defaultdict
from scipy.stats import fisher_exact

def count_activated_mirnas(file_path):
    """Counts the activated miRNAs from the file and returns a dictionary with miRNA names as keys."""
    mirna_counts = defaultdict(int)

    # Open the file and process each line
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line by spaces (or tabs) to extract miRNA name and status
            parts = line.strip().split()
            if len(parts) == 2:
                mirna_name, status = parts
                if status == 'activated':
                    mirna_counts[mirna_name] += 1
    return mirna_counts

def read_lengths(file_path):
    """Reads miRNA lengths from the file and returns a dictionary with miRNA names as keys."""
    mirna_lengths = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Split by comma to extract miRNA name and length
            parts = line.strip().split(', ')
            if len(parts) == 2:
                mirna_name = parts[0]
                length = int(parts[1])
                mirna_lengths[mirna_name] = length
    return mirna_lengths

def merge_mirna_data(activated_counts_file, lengths_file, output_file2):
    """Merges miRNA activated counts and lengths, performs calculations and Fisher's test, then writes the result to a TSV file."""
    # Get the activated counts from the file
    mirna_activated_counts = count_activated_mirnas(activated_counts_file)

    # Get the miRNA lengths from the file
    mirna_lengths = read_lengths(lengths_file)

    # Open the output file in write mode
    with open(output_file2, 'w') as output:
        # Write the header to the output file
        output.write(f"miRNA Name\tcount\tlen\tprocessed_len\terror\ttest\tFisher Result\n")

        # Only keep entries where both activated counts and lengths are available
        common_mirnas = set(mirna_activated_counts.keys()).intersection(set(mirna_lengths.keys()))

        # Loop over the miRNAs and output the information
        for mirna in common_mirnas:
            activated_count = mirna_activated_counts[mirna]
            length = mirna_lengths[mirna]

            # Calculate: Length * 3 - Activated Count
            processed_len = (length * 3) - activated_count

            # Static values for additional columns
            error = 3.86
            test = 189.14

            # Prepare the contingency table for Fisher's test
            contingency_table = [[activated_count, processed_len], [error, test]]

            # Perform Fisher's exact test
            _, p_value = fisher_exact(contingency_table)

            # Determine the result based on the p-value
            fisher_result = "I" if p_value < 0.05 else "IS"

            # Write the final row to the output file in TSV format
            output.write(f"{mirna}\t{activated_count}\t{length}\t{processed_len}\t{error}\t{test}\t{fisher_result}\n")

# Specify the file names
activated_counts_file = user_file + '_temp_activated_hairpins.txt'
lengths_file = user_file2
output_file2 = user_file + '_temp_IS_class.tsv'

# Call the merge function to write the results to a TSV file
merge_mirna_data(activated_counts_file, lengths_file, output_file2)




### ESTRACT DEACTIVATED HAIRPINS

input_file_dr = user_file #"temp_classes_table.tsv"
output_file_dr = user_file + "_temp_deactivated_hairpins.txt"

# Open the input file for reading and output file for writing
with open(input_file_dr, "r") as infile, open(output_file_dr, "w") as outfile:
    # Iterate through each line in the file
    for line in infile:
        # Check if the line contains the desired status
        if "\tdeactivated" in line:
            # Split the line into columns by tab
            columns = line.strip().split("\t")
            # Extract the first and last columns
            first_column = columns[0]
            last_column = columns[-1]
            # Write the extracted columns to the output file
            outfile.write(f"{first_column}\t{last_column}\n")



### DO THE FISHER FOR DR CLASS

def count_deactivated_mirnas(file_path):
    """Counts the deactivated miRNAs from the file and returns a dictionary with miRNA names as keys."""
    mirna_counts = defaultdict(int)
    # Open the file and process each line
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line by spaces (or tabs) to extract miRNA name and status
            parts = line.strip().split()
            if len(parts) == 2:
                mirna_name, status = parts
                if status == 'deactivated':
                    mirna_counts[mirna_name] += 1
    return mirna_counts

def read_lengths(file_path):
    """Reads miRNA lengths from the file and returns a dictionary with miRNA names as keys."""
    mirna_lengths = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Split by comma to extract miRNA name and length
            parts = line.strip().split(', ')
            if len(parts) == 2:
                mirna_name = parts[0]
                length = int(parts[1])
                mirna_lengths[mirna_name] = length
    return mirna_lengths

def merge_mirna_data(deactivated_counts_file, lengths_file, output_file3):
    """Merges miRNA activated counts and lengths, performs calculations and Fisher's test, then writes the result to a TSV file."""
    # Get the activated counts from the file
    mirna_deactivated_counts = count_deactivated_mirnas(deactivated_counts_file)

    # Get the miRNA lengths from the file
    mirna_lengths = read_lengths(lengths_file)

    # Open the output file in write mode
    with open(output_file3, 'w') as output:
        # Write the header to the output file
        output.write(f"miRNA Name\tcount\tlen\tprocessed_len\terror\ttest\tFisher Result\n")

        # Only keep entries where both activated counts and lengths are available
        common_mirnas = set(mirna_deactivated_counts.keys()).intersection(set(mirna_lengths.keys()))

        # Loop over the miRNAs and output the information
        for mirna in common_mirnas:
            deactivated_count = mirna_deactivated_counts[mirna]
            length = mirna_lengths[mirna]

            # Calculate: Length * 3 - Activated Count
            processed_len = (length * 3) - deactivated_count

            # Static values for additional columns
            error = 13.51
            test = 179.49

            # Prepare the contingency table for Fisher's test
            contingency_table = [[deactivated_count, processed_len], [error, test]]

            # Perform Fisher's exact test
            _, p_value = fisher_exact(contingency_table)

            # Determine the result based on the p-value
            fisher_result = "D" if p_value < 0.05 else "DR"

            # Write the final row to the output file in TSV format
            output.write(f"{mirna}\t{deactivated_count}\t{length}\t{processed_len}\t{error}\t{test}\t{fisher_result}\n")

# Specify the file names
deactivated_counts_file = user_file + '_temp_deactivated_hairpins.txt'
lengths_file = user_file2
output_file3 = user_file + '_temp_DR_class.tsv'
# Call the merge function to write the results to a TSV file
merge_mirna_data(deactivated_counts_file, lengths_file, output_file3)
