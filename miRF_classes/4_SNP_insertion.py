import os
import sys

if __name__=="__main__":
   user_file = sys.argv[1]
   user_file2 = sys.argv[2]
#file_path = user_file


# convert the T to U
#input_file = "fasta.fa"
#output_file = "fasta_U.fa"

# Open the input file for reading
#with open(user_file2, 'r') as infile:
    # Open the output file for writing
 #   with open(output_file, 'w') as outfile:
  #      for line in infile:
   #         # If the line starts with ">", it's a header; write it directly to the output file
    #        if line.startswith(">"):
     #           outfile.write(line)
      #      else:
       #         # Convert 'T' to 'U' in the sequence lines
        #        converted_sequence = line.replace('T', 'U')
         #       outfile.write(converted_sequence)


def convert_bases(seq):
    # Define a dictionary to map each base to the other 3 possible bases
    base_dict = {'A': 'CGU', 'C': 'AGU', 'G': 'ACU', 'U': 'ACG'}
    # Initialize an empty list to store the converted sequences
    converted_seqs = []
    # Loop through each character of the sequence
    for i in range(len(seq)):
        # Loop through each possible replacement base
        for j in base_dict[seq[i]]:
            # Replace the base in the sequence with the replacement base
            new_seq = seq[:i] + j + seq[i+1:]
            # Add the converted sequence to the list
            converted_seqs.append(new_seq)
    # Return the list of converted sequences
    return converted_seqs

import glob

# Read the miRNA sequences and their respective headers
mir_dict = {}
#for filename in glob.glob(user_file + "*_final_hairpins.txt"):
with open(user_file2, 'r') as f:
    RNAfold_output = f.readlines()
    header = []
    seq = []
    dots_brackets1 = []
    for i in range(0, len(RNAfold_output), 3):
        header.append(RNAfold_output[i].strip("\n"))
        seq.append(RNAfold_output[i + 1])
        dots_brackets1.append(RNAfold_output[i + 2])
    for e in range(len(header)):
        mir_dict[(header[e].strip("\n"), seq[e].strip("\n"))] = dots_brackets1[e]
f.close()

# Read the complete miRNA sequences from the file
complete_sequences = {}
#for filec in glob.glob(user_file + ".fa"):
with open(user_file, 'r') as f:
    lines = f.readlines()
header = ""
for line in lines:
    line = line.strip()
    if line.startswith(">"):
        header = line
        complete_sequences[header] = ""
    else:
        complete_sequences[header] += line
f.close()

# Convert the bases in each miRNA sequence and write the results to a file
with open(user_file + '_temp_SNPpot_SEQinsterted_fasta.txt', 'w') as f:
    for (header, seq) in mir_dict.keys():
        # Get the complete sequence for the current miRNA
        complete_seq = complete_sequences.get(header, "")
        # Convert the bases in the sequence
        converted_seqs = convert_bases(seq)
        # Write the converted sequences to the output file
        for converted_seq in converted_seqs:
            # in the header of each mirna can be added the type of mutation and its relative position
            # Get the mutation information
            mutation_info = ""
            for i in range(len(seq)):
                if seq[i] != converted_seq[i]:
                    mutation_info += f'{seq[i]}-{converted_seq[i]}, pos={i+1} '
            f.write(header + ' ' + mutation_info.strip() + '\n')
            if complete_seq:
                f.write(complete_seq[:complete_seq.index(seq)] + converted_seq + complete_seq[complete_seq.index(seq) +
                                                                                              len(seq):] + '\n')
            else:
                f.write(converted_seq + '\n')
f.close()

# Run RNAfold command with the final output file as input
import subprocess
output_file = user_file + '_temp_SNPpot_SEQinsterted_fasta.txt'
rnafold_output_file = user_file2 + '_potential_RNAfold.txt'  # Define where RNAfold output will be saved
try:
    # Construct the RNAfold command
    command = ["RNAfold", "-p", "-d2", "--noLP", "--noDP", "--noPS", "--jobs=40", output_file] #command = ["RNAfold", "-p", "-d2", "--noLP", "--noDP", "--no
PS", output_file]
    # Open the RNAfold output file for writing
    with open(rnafold_output_file, "w") as rnafold_out:
        # Run the RNAfold command
        subprocess.run(command, stdout=rnafold_out, check=True)
    # Print success message
    print(f"RNAfold output saved to {rnafold_output_file}")
except subprocess.CalledProcessError as e:
    # Print error message if the subprocess fails
    print(f"RNAfold execution failed: {e}")
