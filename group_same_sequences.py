def read_fasta(file_path):
    sequences = {}
    with open(file_path, "r") as file:
        current_id = None
        current_seq = []

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    seq_str = "".join(current_seq)
                    if seq_str in sequences:
                        sequences[seq_str].append(current_id)
                    else:
                        sequences[seq_str] = [current_id]
                current_id = line
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            seq_str = "".join(current_seq)
            if seq_str in sequences:
                sequences[seq_str].append(current_id)
            else:
                sequences[seq_str] = [current_id]

    return sequences

def write_grouped_fasta(output_path, sequences):
    with open(output_path, "w") as file:
        for seq, headers in sequences.items():
            file.write("|".join(headers) + "\n")
            file.write(seq + "\n")

# Esegui il codice
file_path = "anim_panpan_U.fa" 
output_path = "panpan_grouped.fa" 
sequences = read_fasta(file_path)
write_grouped_fasta(output_path, sequences)
