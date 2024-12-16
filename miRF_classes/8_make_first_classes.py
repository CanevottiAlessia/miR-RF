import sys
import os

if __name__=="__main__":
    user_file = sys.argv[1]

first_input = user_file #"temp_classes_table.tsv"
first_output = user_file + "_col2_classes_temp.txt"

with open(first_input, "r") as f_in, open(first_output, "w") as f_out:
    for line in f_in:
        columns = line.strip().split("\t")
        first_col = columns[0]
        last_col = columns[-1]
        f_out.write(f"{first_col}\t{last_col}\n")



import pandas as pd

input_file = user_file + "_col2_classes_temp.txt"
output_file = user_file + "_temp.txt"

mirna_dict = {}

with open(input_file, "r") as file:
    for line in file:
        if line.startswith(">"):
            parts = line.strip().split("\t")
            mirna_name = parts[0]
            mirna_status = parts[1]
            if mirna_name not in mirna_dict:
                mirna_dict[mirna_name] = set()
            mirna_dict[mirna_name].add(mirna_status)

with open(output_file, "w") as file:
    for mirna_name, statuses in mirna_dict.items():
        if "deactivated" in statuses:
            file.write(f"{mirna_name}\tDispensable\n")
        if "activated" in statuses:
            file.write(f"{mirna_name}\tInducible\n")
        if "noImpact" in statuses:
            file.write(f"{mirna_name}\tSpurious\n")
        else:
            file.write(f"{mirna_name}\tResilient\n")

def classify_mirna(mirbase_output_counts):
    mirna_dict = {}

    # Parse the input and populate the mirna_dict
    for line in mirbase_output_counts:
        if line.strip():  # Use strip to remove leading/trailing whitespaces, if any
            parts = line.split('\t')
            mirna_name = parts[0]
            classification = parts[1].strip()  # Strip to remove whitespaces from the classification

            if mirna_name not in mirna_dict:
                mirna_dict[mirna_name] = set()

            mirna_dict[mirna_name].add(classification)

    # Classify miRNAs based on the criteria
    classified_mirnas = {}

    for mirna_name, classifications in mirna_dict.items():
        if 'Dispensable' in classifications:
            classified_mirnas[mirna_name] = 'Dispensable'
        elif 'Inducible' in classifications:
            classified_mirnas[mirna_name] = 'Inducible'
        else:
            classified_mirnas[mirna_name] = 'Resilient' if 'Resilient' in classifications else 'Spurious'

    return classified_mirnas

# Example usage
mirbase_output_counts_file = open(user_file + "_temp.txt", "r") #("mirgene_output_counts", "r")
mirbase_output_counts = mirbase_output_counts_file.readlines()
mirbase_output_counts_file.close()

result = classify_mirna(mirbase_output_counts)

# for mirna, classification in result.items():
#     print(f"{mirna}: {classification}")

# I want the counts for each class

def classify_mirna(mirbase_output_counts):
    mirna_dict = {}
    class_counts = {'Dispensable': 0, 'Inducible': 0, 'Resilient': 0, 'Spurious': 0}

    # Parse the input and populate the mirna_dict
    for line in mirbase_output_counts:
        if line.strip():  # Use strip to remove leading/trailing whitespaces, if any
            parts = line.split('\t')
            mirna_name = parts[0]
            classification = parts[1].strip()  # Strip to remove whitespaces from the classification

            if mirna_name not in mirna_dict:
                mirna_dict[mirna_name] = set()

            mirna_dict[mirna_name].add(classification)

    # Classify miRNAs based on the criteria and count occurrences
    for mirna_name, classifications in mirna_dict.items():
        if 'Dispensable' in classifications:
            class_counts['Dispensable'] += 1
        elif 'Inducible' in classifications:
            class_counts['Inducible'] += 1
        else:
            class_counts['Resilient'] += 1 if 'Resilient' in classifications else 0
            class_counts['Spurious'] += 1 if 'Spurious' in classifications else 0

    return class_counts

# Example usage
mirbase_output_counts_file = open(user_file + "_temp.txt", "r")
mirbase_output_counts = mirbase_output_counts_file.readlines()
mirbase_output_counts_file.close()

result = classify_mirna(mirbase_output_counts)

#for classification, count in result.items():
    #print(f"{classification}: {count}")



# Define a dictionary to store the classification for each miRNA
mirna_classification = {}

# Read and process the input file
with open(user_file + '_temp.txt', 'r') as file:
    for line in file:
        # Split the line into miRNA name and class
        mirna, mirna_class = line.strip().split()
        # Update the classification based on the rules
        if mirna not in mirna_classification:
            mirna_classification[mirna] = mirna_class
        else:
            # Update classification based on rules
            if mirna_class == 'Inducible':
                mirna_classification[mirna] = 'Inducible'
            elif mirna_class == 'Spurious' and mirna_classification[mirna] != 'Inducible':
                mirna_classification[mirna] = 'Spurious'
            elif mirna_class == 'Dispensable' and mirna_classification[mirna] != 'Inducible' and mirna_classification[mirna] != 'Spurious':
                mirna_classification[mirna] = 'Dispensable'
import os

# Write the results to a file
with open(user_file + '_temp_classes_mirnas.txt', 'w') as output_file:
    for mirna, classification in mirna_classification.items():
        output_file.write(f"{mirna}\t{classification}\n")

#os.remove('temp.txt')
