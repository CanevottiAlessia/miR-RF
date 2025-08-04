import os
import sys
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# === 1. PARSE ARGOMENTI ===
if __name__ == "__main__":
    user_file = sys.argv[1]   # Esempio: "temp_classes_table.tsv"
    user_file2 = sys.argv[2]  # Esempio: "mirna_lengths.txt"

# === 2. ESTRAI HAIRPIN ATTIVATI ===
input_file = user_file
output_file_activated = user_file + "_temp_activated_hairpins.txt"

with open(input_file, "r") as infile, open(output_file_activated, "w") as outfile:
    for line in infile:
        if "\tactivated" in line:
            columns = line.strip().split("\t")
            first_column = columns[0]
            last_column = columns[-1]
            outfile.write(f"{first_column}\t{last_column}\n")

# === 3. ESTRAI HAIRPIN DISATTIVATI ===
output_file_deactivated = user_file + "_temp_deactivated_hairpins.txt"

with open(input_file, "r") as infile, open(output_file_deactivated, "w") as outfile:
    for line in infile:
        if "\tdeactivated" in line:
            columns = line.strip().split("\t")
            first_column = columns[0]
            last_column = columns[-1]
            outfile.write(f"{first_column}\t{last_column}\n")

# === 4. FUNZIONI COMUNI ===
def count_mirnas_by_status(file_path, status_to_match):
    mirna_counts = defaultdict(int)
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                mirna_name, status = parts
                if status == status_to_match:
                    mirna_counts[mirna_name] += 1
    return mirna_counts

def read_lengths(file_path):
    mirna_lengths = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(', ')
            if len(parts) == 2:
                mirna_name = parts[0]
                length = int(parts[1])
                mirna_lengths[mirna_name] = length
    return mirna_lengths

def perform_fisher_analysis(mirna_counts, mirna_lengths, error_val, test_val, output_path, label_positive, label_negative):
    common_mirnas = set(mirna_counts.keys()).intersection(mirna_lengths.keys())
    results = []
    p_values = []

    for mirna in common_mirnas:
        count = mirna_counts[mirna]
        length = mirna_lengths[mirna]
        processed_len = (length * 3) - count
        contingency_table = [[count, processed_len], [error_val, test_val]]
        _, p_value = fisher_exact(contingency_table, alternative="greater")
        results.append((mirna, count, length, processed_len, error_val, test_val, p_value))
        p_values.append(p_value)

    # Calcola FDR con metodo Benjamini-Hochberg
    _, fdr_values, _, _ = multipletests(p_values, method='fdr_bh')

    # Scrivi risultati
    with open(output_path, 'w') as output:
        output.write("miRNA Name\tcount\tlen\tprocessed_len\terror\ttest\tp-value\tFDR\tFisher Result\n")
        for i, (mirna, count, length, proc_len, error, test, pval) in enumerate(results):
            fdr = fdr_values[i]
            result_label = label_positive if fdr < 0.05 else label_negative
            output.write(f"{mirna}\t{count}\t{length}\t{proc_len}\t{error}\t{test}\t{pval:.4g}\t{fdr:.4g}\t{result_label}\n")

# === 5. ESEGUI FISHER PER ATTIVATI ===
activated_counts = count_mirnas_by_status(output_file_activated, "activated")
mirna_lengths = read_lengths(user_file2)
output_file_IS = user_file + "_temp_IS_class.tsv"

perform_fisher_analysis(
    mirna_counts=activated_counts,
    mirna_lengths=mirna_lengths,
    error_val=3.86,
    test_val=189.14,
    output_path=output_file_IS,
    label_positive="I",
    label_negative="IS"
)

# === 6. ESEGUI FISHER PER DISATTIVATI ===
deactivated_counts = count_mirnas_by_status(output_file_deactivated, "deactivated")
output_file_DR = user_file + "_temp_DR_class.tsv"

perform_fisher_analysis(
    mirna_counts=deactivated_counts,
    mirna_lengths=mirna_lengths,
    error_val=13.51,
    test_val=179.49,
    output_path=output_file_DR,
    label_positive="D",
    label_negative="DR"
)
