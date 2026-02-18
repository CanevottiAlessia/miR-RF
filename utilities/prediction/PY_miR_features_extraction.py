import itertools
import pandas as pd
import re
import sys
from pathlib import Path
import string
import random


def generate_random_string(length=6):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for _ in range(length))


if __name__ == "__main__":
    # argv[1] = input file
    # argv[2] = run workdir (temp files directory)
    # argv[3] = final feature table path (tsv) -> used by R
    if len(sys.argv) != 4:
        print("Usage: python3 PY_miR_features_extraction.py <INPUT> <RUN_WORKDIR> <FEATURE_TABLE_OUTPUT>")
        sys.exit(1)

    user_file = sys.argv[1]
    run_workdir = Path(sys.argv[2]).resolve()
    features_table_path = Path(sys.argv[3]).resolve()

    run_workdir.mkdir(parents=True, exist_ok=True)

    safe_base = Path(user_file).name
    file_path = user_file

    # -------------------------
    # STEP 0: preprocess headers (replace tabs in header lines)
    # -------------------------
    output_file = run_workdir / f"header_file_for_{safe_base}"

    with open(file_path, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    for line in lines:
        if line.startswith('>'):
            modified_lines.append(line.replace('\t', ' '))
        else:
            modified_lines.append(line)

    with open(output_file, 'w') as file:
        file.writelines(modified_lines)

    # Read preprocessed file
    with open(output_file, "r") as file_output_rnafold:
        RNAfold_output = file_output_rnafold.readlines()

    if len(RNAfold_output) % 6 != 0:
        print("Error: The input file does not have the expected number of lines per entry.")
        sys.exit(1)

    # -------------------------
    # EXTRACT ENERGIES PARAMETERS
    # -------------------------
    header = []
    seq = []
    dots_brackets1 = []
    dots_brackets2 = []
    dots_brackets3 = []
    energy_MFE = []

    for i in range(0, len(RNAfold_output), 6):
        header.append(RNAfold_output[i].strip("\n"))
        seq.append(RNAfold_output[i + 1])
        dots_brackets1.append(RNAfold_output[i + 2])
        dots_brackets2.append(RNAfold_output[i + 3])
        dots_brackets3.append(RNAfold_output[i + 4])
        energy_MFE.append(RNAfold_output[i + 5])

    mirna_data = []
    for i in range(len(header)):
        mirna = {}
        mirna['header'] = str(header[i])

        energy_match = re.search(r'([-+]\d+\.\d+)', dots_brackets1[i])
        if energy_match:
            mirna['energy_value1'] = float(energy_match.group(1))

        energy_match2 = re.search(r'\[([-+]\d+\.\d+)\]', dots_brackets2[i])
        if energy_match2:
            mirna['energy_value2'] = float(energy_match2.group(1))

        energy_match3 = re.search(r'{\s*([-+]?\d+(\.\d+)?)\s*d\s*=\s*([-+]?\d+(\.\d+)?)\s*}', dots_brackets3[i])
        if energy_match3:
            mirna['energy_value3'] = float(energy_match3.group(1))
            mirna['d_value'] = float(energy_match3.group(3))

        freq_match = re.search(r'ensemble\s([\d.]+)', energy_MFE[i])
        if freq_match:
            mirna['frequency'] = float(freq_match.group(1))

        div_match = re.search(r'diversity\s([\d.]+)', energy_MFE[i])
        if div_match:
            mirna['diversity'] = float(div_match.group(1))

        mirna_data.append(mirna)

    df_energies = pd.DataFrame(mirna_data)
    df_energies.set_index('header', inplace=True)
    df_energies = df_energies.drop_duplicates()

    # -------------------------
    # CLEAN NOTATIONS -> processing file (3 lines per entry)
    # -------------------------
    mir_dict = {}
    for e in range(len(header)):
        mir_dict[(header[e].strip("\n"), seq[e].strip("\n"))] = [
            dots_brackets1[e].strip("\n")[0:len(dots_brackets1[e]) - 10],
            dots_brackets2[e].strip("\n").split("[")[0].replace(" ", ""),
            dots_brackets3[e].strip("\n").split("{")[0].replace(" ", "")
        ]

    resulting_string = run_workdir / (generate_random_string() + "_processing_file.txt")
    with open(resulting_string, "w") as file:
        for d in mir_dict:
            my_name = d[0]
            my_seq = d[1]
            my_not = mir_dict[d][0]

            count_x = 0
            count_y = 0

            for x in my_not:
                if x == ".":
                    count_x += 1
                else:
                    break

            my_not_inv = my_not[::-1]
            for y in my_not_inv:
                if y == ".":
                    count_y += 1
                else:
                    break

            new_not = my_not[count_x:len(my_not) - count_y]
            new_seq = my_seq[count_x:len(my_seq) - count_y]

            file.write(my_name + "\n")
            file.write(new_seq + "\n")
            file.write(new_not + "\n")

    # -------------------------
    # EXTRACTING THE LONGEST HAIRPIN
    # -------------------------
    with open(resulting_string, "r") as f:
        RNAfold_output2 = f.readlines()

    header2 = []
    seq2 = []
    dots2 = []
    mir_dict2 = {}

    for i in range(0, len(RNAfold_output2), 3):
        header2.append(RNAfold_output2[i].strip("\n"))
        seq2.append(RNAfold_output2[i + 1])
        dots2.append(RNAfold_output2[i + 2])

    for e in range(len(header2)):
        mir_dict2[(header2[e].strip("\n"), seq2[e].strip("\n"))] = dots2[e]

    patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
    hairpin_lengths = {}

    for d in mir_dict2:
        my_name = d[0]
        my_seq = d[1]
        my_not = mir_dict2[d]

        r = len(re.findall(patterns, my_not))

        if r == 0:
            if str(re.search(patterns, my_not)) == 'None':
                hairpin_lengths[my_name] = [(my_seq), (my_not)]

        elif r == 1:
            r2 = str(re.search(patterns, my_not))
            indexes = r2.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
            first_hairpin = my_not[:int(indexes[0]) + 1].strip("\n")
            new_seq1 = my_seq[0:len(first_hairpin)]
            second_hairpin = my_not[int(indexes[1]) - 1:].strip("\n")
            new_seq2 = my_seq[len(my_not) - len(second_hairpin) - 1:]
            hairpin_lengths[my_name] = [(new_seq1), (first_hairpin), (new_seq2), (second_hairpin)]

        elif r == 2:
            r2 = str(re.search(patterns, my_not))
            indexes = r2.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
            first_hairpin = my_not[:int(indexes[0]) + 1].strip("\n")
            new_seq1 = my_seq[0:len(first_hairpin)]
            second_hairpin = my_not[int(indexes[1]) - 1:].strip("\n")
            new_seq2 = my_seq[len(my_not) - len(second_hairpin) - 1:]
            sec_pattern = str(re.search(patterns, second_hairpin))
            ind = sec_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
            sec_hairpin = second_hairpin[:int(ind[0]) + 1].strip("\n")
            new_seq_sec = my_seq[len(my_not) - len(second_hairpin) - 1:len(my_not) - len(second_hairpin) - 1 + len(sec_hairpin)]
            third_hairpin = second_hairpin[int(ind[1]) - 1:].strip("\n")
            new_seq3 = my_seq[len(my_not) - len(third_hairpin) - 1:]
            hairpin_lengths[my_name] = [(new_seq1), (first_hairpin), (new_seq_sec), (sec_hairpin), (new_seq3), (third_hairpin)]

        elif r == 3:
            r2 = str(re.search(patterns, my_not))
            indexes = r2.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
            first_hairpin = my_not[:int(indexes[0]) + 1].strip("\n")
            new_seq1 = my_seq[0:len(first_hairpin)]
            second_hairpin = my_not[int(indexes[1]) - 1:].strip("\n")
            new_seq2 = my_seq[len(my_not) - len(second_hairpin) - 1:]
            sec_pattern = str(re.search(patterns, second_hairpin))
            ind = sec_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
            sec_hairpin = second_hairpin[:int(ind[0]) + 1].strip("\n")
            new_seq_sec = my_seq[len(my_not) - len(second_hairpin) - 1:len(my_not) - len(second_hairpin) - 1 + len(sec_hairpin)]
            third_hairpin = second_hairpin[int(ind[1]) - 1:].strip("\n")
            new_seq3 = my_seq[len(my_not) - len(third_hairpin) - 1:]
            third_pattern = str(re.search(patterns, third_hairpin))
            inde = third_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
            th_hairpin = third_hairpin[:int(inde[0]) + 1].strip("\n")
            new_seq_th = my_seq[len(my_not) - len(third_hairpin) - 1:len(my_not) - len(third_hairpin) - 1 + len(th_hairpin)]
            four_hairpin = third_hairpin[int(inde[1]) - 1:].strip("\n")
            new_seq4 = my_seq[len(my_not) - len(four_hairpin) - 1:]
            hairpin_lengths[my_name] = [(new_seq1), (first_hairpin), (new_seq_sec), (sec_hairpin),
                                        (new_seq_th), (th_hairpin), (new_seq4), (four_hairpin)]

        else:
            # fallback (r >= 4): keep whole structure
            hairpin_lengths[my_name] = [(my_seq), (my_not)]

    for key in hairpin_lengths:
        hairpin_lengths[key] = [value.rstrip() for value in hairpin_lengths[key]]

    transformed_dict = {}
    for key, value in hairpin_lengths.items():
        tuples_list = transformed_dict.get(key, [])
        for i in range(0, len(value), 2):
            tuples_list.append((value[i], value[i + 1]))
        transformed_dict[key] = tuples_list

    extracted_dict = {}
    for mirna_name, tuples_list in transformed_dict.items():
        longest_tuple = max(tuples_list, key=lambda t: len(t[1]))
        extracted_dict[mirna_name] = longest_tuple

    random_file_name = run_workdir / (generate_random_string() + "_final_hairpins.txt")
    with open(random_file_name, "w") as file:
        for mirna_name, (sequence, structure) in extracted_dict.items():
            file.write(mirna_name + "\n")
            file.write(sequence + "\n")
            file.write(structure + "\n")

    # -------------------------
    # INFO LOOPS
    # -------------------------
    with open(random_file_name, "r") as f:
        RNAfold_output3 = f.readlines()

    header3 = []
    seq3 = []
    dots3 = []
    for i in range(0, len(RNAfold_output3), 3):
        header3.append(RNAfold_output3[i].strip("\n"))
        seq3.append(RNAfold_output3[i + 1].strip("\n"))
        dots3.append(RNAfold_output3[i + 2].strip("\n"))

    re_3_2 = "(\({2}\.{3,}\){2})+"
    file_info_loops = run_workdir / f"info_loops_for_{safe_base}"

    with open(file_info_loops, "w") as file_info:
        for idx in range(len(header3)):
            my_name = header3[idx].replace("\t", ",")
            my_not = dots3[idx]
            find_loops = str(re.search(re_3_2, my_not))
            if find_loops == 'None':
                file_info.write(str(my_name) + "\t0\t0\t0\t0\n")
            else:
                indexes = find_loops.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
                indexes_updated = [int(indexes[0]) + 1, int(indexes[1]) - 2]
                loop_length = indexes_updated[-1] - indexes_updated[0]
                file_info.write(str(my_name) + "\t1\t" + str(indexes_updated[0]) + "\t" + str(indexes_updated[-1]) + "\t" + str(loop_length) + "\n")

    data = []
    with open(file_info_loops, 'r') as file:
        for line in file:
            parts = line.split('\t')
            data.append({
                'Header': parts[0],
                'loop count': int(parts[1]),
                'start loop': int(parts[2]),
                'end loop': int(parts[3]),
                'loop lenght': int(parts[4].strip())
            })
    df_loops = pd.DataFrame(data)
    df_loops.set_index('Header', inplace=True)

    # -------------------------
    # 32 FEATURES  (FIX: U… A… G… C…)
    # -------------------------
    file_32_features = run_workdir / f"info_32_for_{safe_base}"
    mir_dict = {}
    for i in range(len(header3)):
        mir_dict[(header3[i], seq3[i])] = dots3[i]

    with open(file_32_features, "w") as file_32_info:
        for d in mir_dict:
            my_name = d[0]
            my_seq = d[1]
            my_not = mir_dict[d]
            my_new_not = my_not.replace(")", "(")

            e_dict = {("U", "((("): 0, ("A", "((("): 0, ("G", "((("): 0, ("C", "((("): 0,
                      ("U", "((."): 0, ("A", "((."): 0, ("G", "((."): 0, ("C", "((."): 0,
                      ("U", "(.."): 0, ("A", "(.."): 0, ("G", "(.."): 0, ("C", "(.."): 0,
                      ("U", "(.("): 0, ("A", "(.("): 0, ("G", "(.("): 0, ("C", "(.("): 0,
                      ("U", ".(("): 0, ("A", ".(("): 0, ("G", ".(("): 0, ("C", ".(("): 0,
                      ("U", ".(."): 0, ("A", ".(."): 0, ("G", ".(."): 0, ("C", ".(."): 0,
                      ("U", "..("): 0, ("A", "..("): 0, ("G", "..("): 0, ("C", "..("): 0,
                      ("U", "..."): 0, ("A", "..."): 0, ("G", "..."): 0, ("C", "..."): 0}

            for i in range(0, len(my_seq) - 3):
                tre_seq = my_seq[i:i + 3]
                center_seq = tre_seq[1]
                tre_not = my_new_not[i:i + 3]
                tupla = (center_seq, tre_not)
                if tupla not in e_dict:
                    e_dict[tupla] = 0
                else:
                    e_dict[tupla] += 1

            dict_values = list(e_dict.values())
            tab_dict_values = str(dict_values).replace(",", "\t").replace("[", " ").replace("]", "")
            len_hairpin = len(my_not)
            file_32_info.write(my_name + "\t" + tab_dict_values + "\t" + str(len_hairpin) + "\n")

    data2 = []
    with open(file_32_features, 'r') as file:
        for line in file:
            parts = line.split('\t')
            my_name = parts[0]
            features = list(map(int, parts[1:-1]))
            hairpin_length = int(parts[-1].strip())

            data2.append({
            "Header": my_name,
            "U(((": features[0],
            "A(((": features[1],
            "G(((": features[2],
            "C(((": features[3],
            "U((.": features[4],
            "A((.": features[5],
            "G((.": features[6],
            "C((.": features[7],
            "U(..": features[8],
            "A(..": features[9],
            "G(..": features[10],
            "C(..": features[11],
            "U(.(": features[12],
            "A(.(": features[13],
            "G(.(": features[14],
            "C(.(": features[15],
            "U.((": features[16],
            "A.((": features[17],
            "G.((": features[18],
            "C.((": features[19],
            "U.(.": features[20],
            "A.(.": features[21],
            "G.(.": features[22],
            "C.(.": features[23],
            "U..(": features[24],
            "A..(": features[25],
            "G..(": features[26],
            "C..(": features[27],
            "U…": features[28],
            "A…": features[29],
            "G…": features[30],
            "C…": features[31],
            "hairpin length": hairpin_length
            })

    df_32 = pd.DataFrame(data2)
    df_32.set_index('Header', inplace=True)

    # -------------------------
    # (resto del tuo pipeline)
    # Per brevità: qui sotto rimane identico a prima (nucleotidi, base pairing, bulges, concat, save)
    # -------------------------

    # COUNT NUCLEOTIDES
    file_nucleotides_features = run_workdir / f"info_nucleotides_for_{safe_base}"
    with open(file_nucleotides_features, "w") as file_nucl:
        for i in range(len(header3)):
            my_name = header3[i]
            my_seq = seq3[i]
            nt_dict = {"UU": 0, "AA": 0, "GG": 0, "CC": 0,
                       "UC": 0, "AC": 0, "GC": 0, "CG": 0,
                       "UG": 0, "AG": 0, "GA": 0, "GU": 0,
                       "UA": 0, "AU": 0, "CA": 0, "CU": 0}
            for j in range(0, len(my_seq) - 2):
                di = my_seq[j:j + 2]
                if di in nt_dict:
                    nt_dict[di] += 1
            dict_values = list(nt_dict.values())
            tab_dict_values = str(dict_values).replace(",", "\t").replace("[", "").replace("]", "")
            file_nucl.write(my_name + "\t" + str(tab_dict_values) + "\t0\n")

    data3 = []
    with open(file_nucleotides_features, 'r') as file:
        for line in file:
            parts = line.split('\t')
            header0 = parts[0]
            features = list(map(int, parts[1:-1]))
            real_miRNA = int(parts[-1].strip())
            data3.append({
                'Header': header0,
                'count UU': features[0],
                'count AA': features[1],
                'count GG': features[2],
                'count CC': features[3],
                'count UC': features[4],
                'count AC': features[5],
                'count GC': features[6],
                'count CG': features[7],
                'count UG': features[8],
                'count AG': features[9],
                'count GA': features[10],
                'count GU': features[11],
                'count UA': features[12],
                'count AU': features[13],
                'count CA': features[14],
                'count CU': features[15],
                'real miRNA': real_miRNA,
            })
    df_nucleotides = pd.DataFrame.from_dict(data3)
    df_nucleotides.set_index('Header', inplace=True)

# BASE PAIRING
def load_loop_info(input_file_path):
    loop_info = {}
    with open(input_file_path, 'r') as file:
        for col in file:
            if not col.startswith("#"):
                col = col.strip().split("\t")
                mirna_name = col[0]
                loop_start = int(col[2]) if col[2] != 'NA' else None
                loop_end = int(col[3]) if col[3] != 'NA' else None
                loop_info[mirna_name] = (loop_start, loop_end)
    return loop_info
def count_bases_before_loop_start(header, seq, dots_brackets1, loop_info):
    mirna_counts = {}
    for e in range(len(header)):
        my_name = header[e]
        my_seq = seq[e]
        my_not = dots_brackets1[e]
        loop_start, loop_end = loop_info.get(my_name, (None, None))  # Get loop start and end coordinates
        base_counts = {'A': {'paired': 0, 'unpaired': 0},
                       'U': {'paired': 0, 'unpaired': 0},
                       'G': {'paired': 0, 'unpaired': 0},
                       'C': {'paired': 0, 'unpaired': 0}}
        for j in range(len(my_seq)):
            base = my_seq[j]
            notation = my_not[j]
            if loop_start is not None and loop_end is not None:
                if loop_start <= j <= loop_end:
                    continue  # Skip bases inside the loop
            # Add a check to see if the base is a valid one (A, U, G, C)
            if base not in ['A', 'U', 'G', 'C']:
                print(f"Invalid base '{base}' at position {j} for miRNA {my_name}")
                continue  # Skip this base and continue to the next one
            if notation == '(' or notation == ')':
                base_counts[base]['paired'] += 1
            elif notation == '.':
                base_counts[base]['unpaired'] += 1
        mirna_counts[my_name] = base_counts
    return(mirna_counts)

# Example usage:
RNAfold_out = open(random_file_name, "r")
RNAfold_output = RNAfold_out.readlines()
header = []
seq = []
dots_brackets1 = []
for i in range(0, len(RNAfold_output), 3):
    header.append(RNAfold_output[i].strip("\n"))
    seq.append(RNAfold_output[i + 1].strip("\n"))
    dots_brackets1.append(RNAfold_output[i + 2].strip("\n"))
loop_info = load_loop_info(file_info_loops)
mirna_counts = count_bases_before_loop_start(header, seq, dots_brackets1, loop_info)
df_pairs = pd.DataFrame(columns=['C_paired', 'A_paired', 'G_paired', 'U_paired', 'C_unpaired', 'A_unpaired', 'G_unpaired', 'U_unpaired'])
mirna_rows = []
for mirna, counts in mirna_counts.items():
    row = {
        'mirna_name': mirna,
        'C_paired': counts['C']['paired'],
        'A_paired': counts['A']['paired'],
        'G_paired': counts['G']['paired'],
        'U_paired': counts['U']['paired'],
        'C_unpaired': counts['C']['unpaired'],
        'A_unpaired': counts['A']['unpaired'],
        'G_unpaired': counts['G']['unpaired'],
        'U_unpaired': counts['U']['unpaired']
    }
    mirna_rows.append(row)
df_pairs = pd.DataFrame(mirna_rows, columns=['mirna_name', 'C_paired', 'A_paired', 'G_paired', 'U_paired','C_unpaired', 'A_unpaired', 'G_unpaired', 'U_unpaired'])
df_pairs.set_index('mirna_name', inplace=True)
# Create a df where the rows are the mirna names and the columns correspond to the results of:
# count_CG, count_GC, count_UA, count_AU, count_UG, count_GU, count_AA, count_GG, count_UU, count_CC, frequency_couple_CG,
# frequency_couple_GC, frequency_couple_UA, frequency_couple_AU, frequency_couple_UG, frequency_couple_GU, frequency_couple_AA,
# frequency_couple_GG, frequency_couple_UU, frequency_couple_CC
results = {
    'mirna_name': [],
    'count_CG': [],
    'count_GC': [],
    'count_UA': [],
    'count_AU': [],
    'count_UG': [],
    'count_GU': [],
    'count_AA': [],
    'count_GG': [],
    'count_UU': [],
    'count_CC': [],
    'frequency_couple_CG': [],
    'frequency_couple_GC': [],
    'frequency_couple_UA': [],
    'frequency_couple_AU': [],
    'frequency_couple_UG': [],
    'frequency_couple_GU': [],
    'frequency_couple_AA': [],
    'frequency_couple_GG': [],
    'frequency_couple_UU': [],
    'frequency_couple_CC': []
}
for d in mir_dict:
    my_name = d[0]
    my_seq = d[1]
    my_not = mir_dict[d]#[0]
    total_len = len(my_not)
    initial_index = 0
    final_index = total_len - 1
    couples = []  # list that will contain all the base couples found
    indexes_paired_bases = []  # list that will contain all the indexes of all the couples found
    for e in my_not:
        if my_not[initial_index] == "(" and my_not[final_index] == ")":  # there is a pairing
            couple_of_bases = my_seq[initial_index] + my_seq[final_index]  # extract the first base that corresponds to '(' and the last base that corresponds to ')'
            couples.append(couple_of_bases)
            indexes_to_be_added = []  # sub list that will contain the indexes of the 2 bases that compose the pairing
            indexes_to_be_added.append(initial_index)  # index of the base at 5'
            indexes_to_be_added.append(final_index)  # index of the base at 3'
            indexes_paired_bases.append(indexes_to_be_added)
            initial_index = initial_index + 1
            final_index = final_index - 1
        elif my_not[initial_index] == "." and my_not[final_index] == ")":  # no pairing
             # the initial index increases of 1 to move closer to the centre, while the final index remains fixed
             initial_index = initial_index + 1
        elif my_not[initial_index] == "(" and my_not[final_index] == ".":  # no pairing
            # the initial index remains fixed, while the final index should be decrease of 1 to move closer to the centre
            final_index = final_index - 1
        elif my_not[initial_index] == "." and my_not[final_index] == ".":  # no pairing
            # to move closer to the centre the initial index increases of 1, while the final index decreases of 1
            initial_index = initial_index + 1
            final_index = final_index - 1
    # Now we have to divide all the base couples that have been put in 'couples' list:
    CG_couples = ""  # empty string that will contain all the couples CG
    GC_couples = ""
    UA_couples = ""
    AU_couples = ""
    UG_couples = ""
    GU_couples = ""
    AA_couples = ""
    GG_couples = ""
    UU_couples = ""
    CC_couples = ""
    for couple in couples:  # put every couple of bases found in the corresponding string
        if couple == "CG":
            CG_couples = CG_couples + couple
        elif couple == "GC":
            GC_couples = GC_couples + couple
        elif couple == "UA":
            UA_couples = UA_couples + couple
        elif couple == "AU":
            AU_couples = AU_couples + couple
        elif couple == "UG":
            UG_couples = UG_couples + couple
        elif couple == "GU":
            GU_couples = GU_couples + couple
        elif couple == "AA":
            AA_couples = AA_couples + couple
        elif couple == "GG":
            GG_couples = GG_couples + couple
        elif couple == "UU":
            UU_couples = UU_couples + couple
        elif couple == "CC":
            CC_couples = CC_couples + couple
    count_CG = len(CG_couples) / 2  # count CG occurrences there are in the 'CG_couples' string.
    count_GC = len(GC_couples) / 2
    count_UA = len(UA_couples) / 2
    count_AU = len(AU_couples) / 2
    count_UG = len(UG_couples) / 2
    count_GU = len(GU_couples) / 2
    count_AA = len(AA_couples) / 2
    count_GG = len(GG_couples) / 2
    count_UU = len(UU_couples) / 2
    count_CC = len(CC_couples) / 2
    total_numb_pairings = count_CG + count_GC + count_UA + count_AU + count_UG + count_GU + count_AA + count_GG + count_UU + count_CC  # measure the total number of pairings found
    if total_numb_pairings != 0:
        frequency_couple_CG = count_CG / float(total_numb_pairings)
    else:
        frequency_couple_CG = 0.0
    if total_numb_pairings != 0:
        frequency_couple_GC = count_GC / float(total_numb_pairings)
    else:
        frequency_couple_GC = 0.0
    if total_numb_pairings != 0:
        frequency_couple_UA = count_UA / float(total_numb_pairings)
    else:
        frequency_couple_UA = 0.0
    if total_numb_pairings != 0:
        frequency_couple_AU = count_AU / float(total_numb_pairings)
    else:
        frequency_couple_AU = 0.0
    if total_numb_pairings != 0:
        frequency_couple_UG = count_UG / float(total_numb_pairings)
    else:
        frequency_couple_UG = 0.0
    if total_numb_pairings != 0:
        frequency_couple_GU = count_GU / float(total_numb_pairings)
    else:
        frequency_couple_GU = 0.0
    if total_numb_pairings != 0:
        frequency_couple_AA = count_AA / float(total_numb_pairings)
    else:
        frequency_couple_AA = 0.0
    if total_numb_pairings != 0:
        frequency_couple_GG = count_GG / float(total_numb_pairings)
    else:
        frequency_couple_GG = 0.0
    if total_numb_pairings != 0:
        frequency_couple_UU = count_UU / float(total_numb_pairings)
    else:
        frequency_couple_UU = 0.0
    if total_numb_pairings != 0:
        frequency_couple_CC = count_CC / float(total_numb_pairings)
    else:
        frequency_couple_CC = 0.0
    results['mirna_name'].append(my_name)
    results['count_CG'].append(count_CG)
    results['count_GC'].append(count_GC)
    results['count_UA'].append(count_UA)
    results['count_AU'].append(count_AU)
    results['count_UG'].append(count_UG)
    results['count_GU'].append(count_GU)
    results['count_AA'].append(count_AA)
    results['count_GG'].append(count_GG)
    results['count_UU'].append(count_UU)
    results['count_CC'].append(count_CC)
    results['frequency_couple_CG'].append(frequency_couple_CG)
    results['frequency_couple_GC'].append(frequency_couple_GC)
    results['frequency_couple_UA'].append(frequency_couple_UA)
    results['frequency_couple_AU'].append(frequency_couple_AU)
    results['frequency_couple_UG'].append(frequency_couple_UG)
    results['frequency_couple_GU'].append(frequency_couple_GU)
    results['frequency_couple_AA'].append(frequency_couple_AA)
    results['frequency_couple_GG'].append(frequency_couple_GG)
    results['frequency_couple_UU'].append(frequency_couple_UU)
    results['frequency_couple_CC'].append(frequency_couple_CC)
df_b = pd.DataFrame(results)
df_b.set_index('mirna_name', inplace=True)
df_base_features = pd.concat([df_pairs, df_b], axis=1)


# INFO BULGES
RNAfold_out = open(random_file_name, "r")
RNAfold_output = RNAfold_out.readlines()
header = []
seq = []
dots_brackets1 = []
mir_dict = {}
for i in range(0, len(RNAfold_output), 3):
    header.append(RNAfold_output[i].strip("\n"))
    seq.append(RNAfold_output[i + 1].strip("\n"))
    dots_brackets1.append(RNAfold_output[i + 2].strip("\n"))
for e in range(len(header)):
    mir_dict[(header[e].strip("\n"), seq[e].strip("\n"))] = dots_brackets1[e]
pattern_counts = {}  # Dictionary to store pattern counts
miRNA_order = []  # Maintain the order of miRNAs
bulge_5_prime_counts = {}  # Dictionary to store bulge counts at 5' prime
bulge_3_prime_counts = {}  # Dictionary to store bulge counts at 3' prime
for d in mir_dict:
    my_name = d[0]
    my_seq = d[1]
    my_not = mir_dict[d]
    total_length_seq = len(my_not)  # not useful here
    regex_putative_bulge_5_prime = "\(\.{1,}\("
    regex_putative_bulge_3_prime = "\)\.{1,}\)"
    regex_putative_bulge = regex_putative_bulge_5_prime + "|" + regex_putative_bulge_3_prime
    bul = re.findall(regex_putative_bulge, my_not)
    central_loop = re.search(r"(\({2}\.{3,}\){2})+", my_not)  # not useful here
    indices = []  # List to store the start and end indices of matched patterns   
    for b in bul:
        if b not in pattern_counts:
            pattern_counts[b] = {}  # Create an empty dictionary for the pattern if it's not already present
        if my_name not in pattern_counts[b]:
            pattern_counts[b][my_name] = 1  # Initialize the count for the miRNA in the pattern dictionary
        else:
            pattern_counts[b][my_name] += 1  # Increment the count for the miRNA in the pattern dictionary           
        # Find the first occurrence of the pattern in the RNA secondary structure string
        match = re.search(re.escape(b), my_not)
        # Get the start and end indices of the matched pattern
        start = match.start()
        end = match.end()
        # Store the start and end indices as a tuple in the indices list
        indices.append((start, end))       
        if '(' in b:
            if my_name not in bulge_5_prime_counts:
                bulge_5_prime_counts[my_name] = 1
            else:
                bulge_5_prime_counts[my_name] += 1
        elif ')' in b:
            if my_name not in bulge_3_prime_counts:
                bulge_3_prime_counts[my_name] = 1
            else:
                bulge_3_prime_counts[my_name] += 1                
    if my_name not in miRNA_order:
        miRNA_order.append(my_name)
# Define the desired patterns
patterns = [
    '(.(',
    '(............(',
    '(...(',
    ')....)',
    ').)',
    ')..)',
    ').....)',
    ')......)',
    '(..(',
    ')...)',
    '(.....(',
    '(......(',
    '(....(',
    ').......)',
    '(........(',
    ')..........)',
    ')........)',
    ').........)',
    '(.........(',
    '(.......(',
    '(..........(',
    ').............)',
    ')...............)',
    '(...........(',
    '(.................(',
    ')............)',
    ')...........)',
    ')..................)',
    '(..............(',
    '(...............(',
    '(...................(',
    '(.............(',
    '(..................(',
    '(.....................(',
    ')................)'
]
# Initialize pattern counts to 0
for pattern in patterns:
    pattern_counts.setdefault(pattern, {})
df = pd.DataFrame(pattern_counts)
df['Bulges_5_prime'] = df.index.map(bulge_5_prime_counts)
df['Bulges_3_prime'] = df.index.map(bulge_3_prime_counts)
df = df.reindex(miRNA_order)
df = df.fillna(0)
column_order = ['(.(', '(............(', '(...(', ')....)', ').)', ')..)', ').....)', ')......)',
                '(..(', ')...)', '(.....(', '(......(', '(....(', ').......)', '(........(',
                ')..........)', ')........)', ').........)', '(.........(', '(.......(',
                '(..........(', ').............)', ')...............)', '(...........(',
                '(.................(', ')............)', ')...........)', ')..................)',
                '(..............(', '(...............(', '(...................(', '(.............(',
                '(..................(', '(.....................(', ')................)',
                'Bulges_5_prime', 'Bulges_3_prime']
df_bulges = df[column_order]
# df_bulges.to_csv('30_info_bulges.tsv', sep="\t")

DF_TOTAL = pd.DataFrame(pd.concat([df_loops, df_32, df_nucleotides, df_base_features, df_bulges, df_energies], axis=1))
#DF_TOTAL = pd.concat([df_loops, df_32, df_nucleotides, df_base_features, df_bulges, df_energies], axis=1, join='inner')
DF_TOTAL.to_csv(features_table_path, sep="\t")


import os

# cleanup intermedi (lascia features_table_path)
for f in [output_file, resulting_string, random_file_name, file_info_loops, file_32_features, file_nucleotides_features]:
    try:
        Path(f).unlink()
    except FileNotFoundError:
        pass
