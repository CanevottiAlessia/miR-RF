import itertools
from typing import TextIO
import pandas as pd
import re
import os
import sys


if __name__=="__main__":
    user_file = sys.argv[1]
file_path = user_file
output_file = f"temp_header_file_for_{user_file}"


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


try:
    with open(output_file, "r") as file_output_rnafold:
        RNAfold_output = file_output_rnafold.readlines()
except FileNotFoundError:
    print(f"Error: The input file '{user_file}' was not found.")
    quit()
except IOError as e:
    print(f"An error occurred while opening the file: {e}")
    quit()
if len(RNAfold_output) % 6 != 0:
    print("Error: The input file does not have the expected number of lines per entry.")
    quit()



# EXTRACT ENERGIES PARAMETERS
header = []  # This is for the headers
seq = []  # This is for the sequences
dots_brackets1 = []  # This is for the first dots-brackets notation
dots_brackets2 = []  # This is for the second dots-brackets notation
dots_brackets3 = []  # This is for the third dots-brackets notation
energy_MFE = []  # This is for the energy-MFE values
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
    mirna['header'] = str(header[i])#.replace("\t", ",")).split(",")
    #mirna['header'] = mirna['header'][0]
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
#df_energies['header'] = df_energies['header'].apply(lambda x: '"' + x + '"')
df_energies.set_index('header', inplace=True)
df_energies = df_energies.drop_duplicates()


# I create a function that removes the dots only at the beginning and at the end of each dots_brackets1 notation,
# creating a file that only has the clean notations (clean_seq.out) which will be the new input file
file_output_rnafold = open(output_file, "r")
RNAfold_output = file_output_rnafold.readlines()
header = []  # This is for the headers
seq = []  # This is for the sequences
dots_brackets1 = []  # This is for the first dots-brackets notation
dots_brackets2 = []  # This is for the second dots-brackets notation
dots_brackets3 = []  # This is for the third dots-brackets notation
energy_MFE = []  # This is for the energy-MFE values
mir_dict = {}  # empty dictionary
for i in range(0, len(RNAfold_output), 6):
    header.append(RNAfold_output[i].strip("\n"))
    seq.append(RNAfold_output[i + 1])
    dots_brackets1.append(RNAfold_output[i + 2])
    dots_brackets2.append(RNAfold_output[i + 3])
    dots_brackets3.append(RNAfold_output[i + 4])
    energy_MFE.append(RNAfold_output[i + 5])
for e in range(len(header)):
    mir_dict[(header[e].strip("\n"), seq[e].strip("\n"))] = [
        dots_brackets1[e].strip("\n")[0:len(dots_brackets1[e]) - 10],
        dots_brackets2[e].strip("\n").split("[")[0].replace(" ", ""),
        dots_brackets3[e].strip("\n").split("{")[0].replace(" ", "")
    ]
import string
import random
def generate_random_string(length=6):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for _ in range(length))
resulting_string = generate_random_string() + "_processing_file.txt"
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
        #new_name = str(my_name.replace("\t", ",")).split(",")
        #name = new_name[0]
        #name = str('"' + my_name + '"')
        # Write the data to the file
        file.write(my_name + "\n")
        file.write(new_seq + "\n")
        file.write(new_not + "\n")



# EXTRACTING THE LONGEST HAIRPIN
import re
file_output_rnafold = open(resulting_string, "r")
RNAfold_output = file_output_rnafold.readlines()
header = []
seq = []
dots_brackets1 = []
mir_dict = {}
for i in range(0, len(RNAfold_output), 3):
    header.append(RNAfold_output[i].strip("\n"))
    seq.append(RNAfold_output[i + 1])
    dots_brackets1.append(RNAfold_output[i + 2])
for e in range(len(header)):
    mir_dict[(header[e].strip("\n"), seq[e].strip("\n"))] = dots_brackets1[e]
# file = StringIO()
hairpin_lengths = {}
for d in mir_dict:
    my_name = d[0]
    my_seq = d[1]
    my_not = mir_dict[d]
    patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
    r = len(re.findall(patterns, my_not))
    if r == 0:  # 1 hairpin
        patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
        r = str(re.search(patterns, my_not))
        if r == 'None':
            my_not = my_not
            # print(my_name + "\n" + my_seq + "\n" + my_not.strip("\n"))
            hairpin_lengths[my_name] = [(my_seq), (my_not)]
    if r == 1:  # 2 hairpins
        patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
        r = str(re.search(patterns, my_not))
        indexes = r.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        list_of_first_int = []
        first_index = [indexes[0], ]
        for e in first_index:
            list_of_first_int.append(int(e))
        first_hairpin = my_not[:list_of_first_int[0] + 1].strip("\n")
        new_seq1 = my_seq[0:len(first_hairpin)]  # sequence for the first hairpin (seems OK)
        list_of_last_int = []
        last_index = [indexes[1], ]
        for i in last_index:
            list_of_last_int.append(int(i))
        second_hairpin = my_not[list_of_last_int[0] - 1:].strip("\n")
        new_seq2 = my_seq[
                   len(my_not) - len(second_hairpin) - 1:]  # sequence for the second hairpin (second or sec+third)
        sec_pattern = str(re.search(patterns, second_hairpin))
        hairpin_lengths[my_name] = [(new_seq1), (first_hairpin), (new_seq2), (second_hairpin)]
    if r == 2:  # 3 hairpins
        patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
        r = str(re.search(patterns, my_not))
        indexes = r.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        list_of_first_int = []
        first_index = [indexes[0], ]
        for e in first_index:
            list_of_first_int.append(int(e))
        first_hairpin = my_not[:list_of_first_int[0] + 1].strip("\n")
        new_seq1 = my_seq[0:len(first_hairpin)]  # sequence for the first hairpin
        list_of_last_int = []
        last_index = [indexes[1], ]
        for i in last_index:
            list_of_last_int.append(int(i))
        second_hairpin = my_not[list_of_last_int[0] - 1:].strip("\n")
        new_seq2 = my_seq[
                   len(my_not) - len(second_hairpin) - 1:]  # sequence for the second hairpin (second or sec+third)
        sec_pattern = str(re.search(patterns, second_hairpin))
        ind = sec_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        l_of_first_ind = []
        f_ind = [ind[0], ]
        for c in f_ind:
            l_of_first_ind.append(int(c))
        sec_hairpin = second_hairpin[:l_of_first_ind[0] + 1].strip("\n")
        new_seq_sec = my_seq[len(my_not) - len(second_hairpin) - 1:len(my_not) - len(second_hairpin) - 1 + len(
            sec_hairpin)]
        l_of_last_ind = []
        l_ind = [ind[1], ]
        for t in l_ind:
            l_of_last_ind.append(int(t))
        third_hairpin = second_hairpin[l_of_last_ind[0] - 1:].strip("\n")
        new_seq3 = my_seq[len(my_not) - len(third_hairpin) - 1:]
        hairpin_lengths[my_name] = [(new_seq1), (first_hairpin), (new_seq_sec), (sec_hairpin), (new_seq3),
                                    (third_hairpin)]
    if r == 3:  # 4 hairpins
        patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
        r = str(re.search(patterns, my_not))
        indexes = r.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        list_of_first_int = []
        first_index = [indexes[0], ]
        for e in first_index:
            list_of_first_int.append(int(e))
        first_hairpin = my_not[:list_of_first_int[0] + 1].strip("\n")
        new_seq1 = my_seq[0:len(first_hairpin)]  # sequence for the first hairpin
        list_of_last_int = []
        last_index = [indexes[1], ]
        for i in last_index:
            list_of_last_int.append(int(i))
        second_hairpin = my_not[list_of_last_int[0] - 1:].strip("\n")
        new_seq2 = my_seq[
                   len(my_not) - len(second_hairpin) - 1:]  # sequence for the second hairpin (second or sec+third)
        sec_pattern = str(re.search(patterns, second_hairpin))
        ind = sec_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        l_of_first_ind = []
        f_ind = [ind[0], ]
        for c in f_ind:
            l_of_first_ind.append(int(c))
        sec_hairpin = second_hairpin[:l_of_first_ind[0] + 1].strip("\n")  # OK
        new_seq_sec = my_seq[len(my_not) - len(second_hairpin) - 1:len(my_not) - len(second_hairpin) - 1 + len(
            sec_hairpin)]
        l_of_last_ind = []
        l_ind = [ind[1], ]
        for t in l_ind:
            l_of_last_ind.append(int(t))
        third_hairpin = second_hairpin[l_of_last_ind[0] - 1:].strip("\n")
        new_seq3 = my_seq[len(my_not) - len(third_hairpin) - 1:]
        third_pattern = str(re.search(patterns, third_hairpin))
        inde = third_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")  # ind =
        l_of_f_ind = []  # l_of_first_ind = []
        fi_ind = [inde[0], ]  # f_ind = [ind[0], ]
        for c in fi_ind:
            l_of_f_ind.append(int(c))
        th_hairpin = third_hairpin[:l_of_f_ind[0] + 1].strip("\n")
        new_seq_th = my_seq[
                     len(my_not) - len(third_hairpin) - 1:len(my_not) - len(third_hairpin) - 1 + len(th_hairpin)]
        l_of_las_ind = []  # l_of_last_ind = []
        li_ind = [inde[1], ]  # l_ind = [ind[1], ]
        for t in li_ind:  # l_ind
            l_of_las_ind.append(int(t))  # l_of_last_ind.append(int(t))
        four_hairpin = third_hairpin[l_of_las_ind[0] - 1:].strip("\n")
        new_seq4 = my_seq[len(my_not) - len(four_hairpin) - 1:]
        hairpin_lengths[my_name] = [(new_seq1), (first_hairpin), (new_seq_sec), (sec_hairpin), (new_seq_th),
                                    (th_hairpin), (new_seq4), (four_hairpin)]
    if r == 4:  # 5 hairpins (6 sequences, 3 miRNAs)
        patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
        r = str(re.search(patterns, my_not))
        indexes = r.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        list_of_first_int = []
        first_index = [indexes[0], ]
        for e in first_index:
            list_of_first_int.append(int(e))
        first_hairpin = my_not[:list_of_first_int[0] + 1].strip("\n")
        new_seq1 = my_seq[0:len(first_hairpin)]  # sequence for the first hairpin
        list_of_last_int = []
        last_index = [indexes[1], ]
        for i in last_index:
            list_of_last_int.append(int(i))
        second_hairpin = my_not[list_of_last_int[0] - 1:].strip("\n")
        new_seq2 = my_seq[len(my_not) - len(second_hairpin) - 1:]
        sec_pattern = str(re.search(patterns, second_hairpin))
        ind = sec_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        l_of_first_ind = []
        f_ind = [ind[0], ]
        for c in f_ind:
            l_of_first_ind.append(int(c))
        sec_hairpin = second_hairpin[:l_of_first_ind[0] + 1].strip("\n")  # OK
        new_seq_sec = my_seq[len(my_not) - len(second_hairpin) - 1:len(my_not) - len(second_hairpin) - 1 + len(
            sec_hairpin)]
        l_of_last_ind = []
        l_ind = [ind[1], ]
        for t in l_ind:
            l_of_last_ind.append(int(t))
        third_hairpin = second_hairpin[l_of_last_ind[0] - 1:].strip("\n")
        new_seq3 = my_seq[len(my_not) - len(third_hairpin) - 1:]
        third_pattern = str(re.search(patterns, third_hairpin))
        inde = third_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")  # ind =
        l_of_f_ind = []
        fi_ind = [inde[0], ]
        for c in fi_ind:
            l_of_f_ind.append(int(c))
        th_hairpin = third_hairpin[:l_of_f_ind[0] + 1].strip("\n")
        new_seq_th = my_seq[
                     len(my_not) - len(third_hairpin) - 1:len(my_not) - len(third_hairpin) - 1 + len(th_hairpin)]
        l_of_las_ind = []
        li_ind = [inde[1], ]
        for t in li_ind:
            l_of_las_ind.append(int(t))
        four_hairpin = third_hairpin[l_of_las_ind[0] - 1:].strip("\n")
        new_seq4 = my_seq[len(my_not) - len(four_hairpin) - 1:]
        four_pattern = str(re.search(patterns, third_hairpin))
        indice = four_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        lista_of_first_ind = []
        fir_ind = [indice[0], ]
        for c in fir_ind:
            lista_of_first_ind.append(int(c))
        fr_hairpin = third_hairpin[:lista_of_first_ind[0] + 1].strip("\n")  # effective third
        new_seq_fr = my_seq[len(my_not) - len(third_hairpin) - 1:len(my_not) - len(third_hairpin) - 1 + len(
            fr_hairpin)]  # OK
        lista_of_last_ind = []
        lis_ind = [indice[1], ]
        for t in lis_ind:
            lista_of_last_ind.append(int(t))
        four_five_hairpin = third_hairpin[lista_of_last_ind[0] - 1:].strip("\n")  # four + five hairpins
        five_pattern = str(re.search(patterns, four_five_hairpin))
        ultimo_index = five_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        lf = []
        fi = [ultimo_index[0], ]
        for x in fi:
            lf.append(int(x))
        eff_four_hairpin = four_five_hairpin[:lf[0] + 1].strip("\n")
        new_seq_eff_four = my_seq[len(my_not) - len(four_hairpin) - 1: len(my_not) - len(four_hairpin) - 1 + len(
            eff_four_hairpin)]
        ll = []
        li = [ultimo_index[1], ]
        for z in li:
            ll.append(int(z))
        eff_five_hairpin = four_five_hairpin[ll[0] - 1:].strip("\n")
        new_seq_eff_five = my_seq[len(my_not) - len(eff_five_hairpin) - 1:]
        hairpin_lengths[my_name] = [(new_seq1), (first_hairpin), (new_seq_sec), (sec_hairpin),
                                    (new_seq_fr), (fr_hairpin), (new_seq_eff_four), (eff_four_hairpin),
                                    (new_seq_eff_five), (eff_five_hairpin)]
    if r == 5:
        patterns = "(\)\()+" "|" "(\){1}\.{1,}\({1})+"
        r = str(re.search(patterns, my_not))
        indexes = r.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        list_of_first_int = []
        first_index = [indexes[0], ]
        for e in first_index:
            list_of_first_int.append(int(e))
        first_hairpin = my_not[:list_of_first_int[0] + 1].strip("\n")
        new_seq1 = my_seq[0:len(first_hairpin)]  # sequence for the first hairpin
        list_of_last_int = []
        last_index = [indexes[1], ]
        for i in last_index:
            list_of_last_int.append(int(i))
        second_hairpin = my_not[list_of_last_int[0] - 1:].strip("\n")
        new_seq2 = my_seq[len(my_not) - len(second_hairpin) - 1:]
        sec_pattern = str(re.search(patterns, second_hairpin))
        ind = sec_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        l_of_first_ind = []
        f_ind = [ind[0], ]
        for c in f_ind:
            l_of_first_ind.append(int(c))
        sec_hairpin = second_hairpin[:l_of_first_ind[0] + 1].strip("\n")  # OK
        new_seq_sec = my_seq[len(my_not) - len(second_hairpin) - 1:len(my_not) - len(second_hairpin) - 1 + len(
            sec_hairpin)]
        l_of_last_ind = []
        l_ind = [ind[1], ]
        for t in l_ind:
            l_of_last_ind.append(int(t))
        third_hairpin = second_hairpin[l_of_last_ind[0] - 1:].strip("\n")
        new_seq3 = my_seq[len(my_not) - len(third_hairpin) - 1:]

        third_pattern = str(re.search(patterns, third_hairpin))
        inde = third_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")  # ind =
        l_of_f_ind = []
        fi_ind = [inde[0], ]
        for c in fi_ind:
            l_of_f_ind.append(int(c))
        th_hairpin = third_hairpin[:l_of_f_ind[0] + 1].strip("\n")
        new_seq_th = my_seq[
                     len(my_not) - len(third_hairpin) - 1:len(my_not) - len(third_hairpin) - 1 + len(th_hairpin)]
        l_of_las_ind = []
        li_ind = [inde[1], ]
        for t in li_ind:
            l_of_las_ind.append(int(t))
        four_hairpin = third_hairpin[l_of_las_ind[0] - 1:].strip("\n")
        new_seq4 = my_seq[len(my_not) - len(four_hairpin) - 1:]
        four_pattern = str(re.search(patterns, third_hairpin))
        indice = four_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        lista_of_first_ind = []
        fir_ind = [indice[0], ]
        for c in fir_ind:
            lista_of_first_ind.append(int(c))
        fr_hairpin = third_hairpin[:lista_of_first_ind[0] + 1].strip("\n")  # effective third
        new_seq_fr = my_seq[len(my_not) - len(third_hairpin) - 1:len(my_not) - len(third_hairpin) - 1 + len(
            fr_hairpin)]  # OK
        lista_of_last_ind = []
        lis_ind = [indice[1], ]
        for t in lis_ind:
            lista_of_last_ind.append(int(t))
        four_five_hairpin = third_hairpin[lista_of_last_ind[0] - 1:].strip("\n")  # four + five + six hairpins
        five_pattern = str(re.search(patterns, four_five_hairpin))
        ultimo_index = five_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        lf = []
        fi = [ultimo_index[0], ]
        for x in fi:
            lf.append(int(x))
        eff_four_hairpin = four_five_hairpin[:lf[0] + 1].strip("\n")  # ok è la 4
        new_seq_eff_four = my_seq[len(my_not) - len(four_hairpin) - 1: len(my_not) - len(four_hairpin) - 1 + len(
            eff_four_hairpin)]
        ll = []
        li = [ultimo_index[1], ]
        for z in li:
            ll.append(int(z))
        eff_five_hairpin = four_five_hairpin[ll[0] - 1:].strip("\n")  # five + six
        new_seq_eff_five = my_seq[len(my_not) - len(eff_five_hairpin) - 1:]  # five + six
        six_pattern = str(re.search(patterns, eff_five_hairpin))
        ultimo_index = six_pattern.split("span=")[1].split(", m")[0].replace("(", "").replace(")", "").split(",")
        lf = []
        fi = [ultimo_index[0], ]
        for x in fi:
            lf.append(int(x))
        effective_five_hairpin = four_five_hairpin[:lf[0] + 1].strip("\n")  # five hairpin ok
        new_seq_effective_five = my_seq[len(my_not) - len(eff_five_hairpin) - 1: len(my_not) - len(
            eff_five_hairpin) - 1 + len(effective_five_hairpin)]
        ll = []
        li = [ultimo_index[1], ]
        for z in li:
            ll.append(int(z))
        six_hairpin = eff_five_hairpin[ll[0] - 1:].strip("\n")
        new_seq_six = my_seq[len(my_not) - len(six_hairpin) - 1:]
        hairpin_lengths[my_name] = [(new_seq1), (first_hairpin),
                                    (new_seq_sec), (sec_hairpin),
                                    (new_seq_fr), (fr_hairpin),
                                    (new_seq_eff_four), (eff_four_hairpin),
                                    (new_seq_effective_five), (effective_five_hairpin),
                                    (new_seq_six), (six_hairpin)]
# Remove '\n' from dictionary values
for key in hairpin_lengths:
    hairpin_lengths[key] = [value.rstrip() for value in hairpin_lengths[key]]
# file.write(str(hairpin_lengths))
transformed_dict = {}
for key, value in hairpin_lengths.items():
    mirna_name = key
    tuples_list = transformed_dict.get(mirna_name, [])
    sequence = ""
    structure = ""
    for i in range(0, len(value), 2):
        sequence = value[i]
        structure = value[i + 1]
        tuples_list.append((sequence, structure))
    transformed_dict[mirna_name] = tuples_list
# Step 4: Extract the tuple with the longest element for each list
#extracted_dict = {}
#for mirna_name, tuples_list in transformed_dict.items():
#    longest_tuple = max(tuples_list, key=lambda t: len(t[1]))
#    extracted_dict[mirna_name] = longest_tuple
#    print(f"{mirna_name}, {len(longest_tuple[1])}")


# Step 4: Extract the tuple with the longest element for each list
extracted_dict = {}
output_file_path = user_file + "_temp_len_hairpins.txt"  # Specify the file name or path

# Open the file for writing
with open(output_file_path, "w") as file:
    for mirna_name, tuples_list in transformed_dict.items():
        longest_tuple = max(tuples_list, key=lambda t: len(t[1]))
        extracted_dict[mirna_name] = longest_tuple
        # Write the output to the file
        print(f"{mirna_name}, {len(longest_tuple[1])}", file=file)

#print(f"Output has been saved to {output_file_path}")




import os
import glob

# Initialize the list to hold the filenames
intermediate_file = []
# Find all files ending with '_final_hairpins.txt' and add them to the list
#  intermediate_file.extend(glob.glob("*_final_hairpins.txt"))
# Find all files ending with '_processing_file.txt' and add them to the list
intermediate_file.extend(glob.glob("*_processing_file.txt"))
intermediate_file.extend(glob.glob("temp_header_file_for_*"))

# Loop through the files and attempt to remove them
for file in intermediate_file:
    try:
        os.remove(file)  # Remove the file
    except:
        pass
