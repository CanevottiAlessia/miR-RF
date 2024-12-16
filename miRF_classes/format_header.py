import os
import sys

if __name__=="__main__":
   user_file = sys.argv[1]
   user_file2 = sys.argv[2]

def replace_spaces_in_headers(file_path, output_path):
    with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Replace spaces and '\t' with '_' in header lines
                line = line.replace(' ', '_').replace('\t', '_')
            outfile.write(line)

# Use
input_file = user_file
output_file = user_file2 
replace_spaces_in_headers(input_file, output_file)
