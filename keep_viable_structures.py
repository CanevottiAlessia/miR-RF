import re

def count_hairpins(structure):
    # Identifies loops
    loops = re.findall(r"\([^()]+\)", structure)
    return len(loops)

def parse_file(input_filename, output_filename):
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    hairpin_exceeding = []
    filtered_lines = []

    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            header = lines[i].strip()
            sequence = lines[i+1].strip()
            structure = lines[i+2].split()[0] # primary structure

            num_hairpins = count_hairpins(structure)
            if num_hairpins > 5:
                hairpin_exceeding.append((header, num_hairpins, structure))
                i += 6  # Goes to the next block
                continue

        filtered_lines.append(lines[i])
        i += 1

    with open(output_filename, 'w') as file:
        file.writelines(filtered_lines)

    return hairpin_exceeding

# Parsing and writing in a new file
input_filename = "RNAfold_output_ok.txt" #"rnafold_PROVA.txt"  # Nome del file originale
output_filename = "RNAfold_output_ok_filt.txt"  # Nome del file con i dati mantenuti
result = parse_file(input_filename, output_filename)

# Stampa i risultati
for header, num_hairpins, structure in result:
    print(f"{header}: {num_hairpins} hairpin trovati e rimossi")
