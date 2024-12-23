import sys
import os
import argparse

def remove_temp_files():
    """Removes all files with 'temp' or 'final_hairpins.' in current directory."""
    for file_name in os.listdir("."):
        if "temp" in file_name or "final_hairpins." in file_name:
            try:
                os.remove(file_name)
                print(f"Removed file: {file_name}")
            except Exception as e:
                print(f"Error removing file {file_name}: {e}")

def runner(rnafold, fasta, base_pred, out):
    try:
        os.system("python3 get_fasta_for_insertion.py " + rnafold)
    except:
        print("Error: Failed to execute. Please, enter the input file name.")
        pass
    intermediate_file = rnafold + "_final_hairpins.txt"
    try:
        os.system("python3 SNP_insertion.py " + fasta + " " + intermediate_file)
    except:
        print("Error: Failed to execute python3 SNP_insertion.py")
        pass
    intermediate_file2 = intermediate_file + "_potential_RNAfold.txt"
    try:
        os.system("python3 2_miR_application.py " + intermediate_file2)
    except:
        print("Error: Failed to execute python3 2_miR_application.py")
        pass
    intermediate_file3 = "features_table_for_" + intermediate_file2 + "_temp_pred_potential.tsv"
    try:
        os.system("python3 merge_table.py " + base_pred + " " + intermediate_file3)
    except:
        print("Error: Failed to execute python3 merge_table.py")
        pass
    intermediate_file4 = base_pred + "_temp_classes_table.tsv"
    try:
        os.system("python3 get_lens.py " + rnafold)
    except OSError as e:
        print("Error: Failed to execute python3 get_lens.py")
        pass
    intermediate_file5 = rnafold + "_temp_len_hairpins.txt"
    try:
        os.system("python3 make_fisher_test.py " + intermediate_file4 + " " + intermediate_file5)
    except:
        print("Error: Failed to execute python3 make_fisher_test.py")
        pass
    try:
        os.system("python3 make_first_classes.py " + intermediate_file4)
    except:
        print("Error: Failed to execute python3 make_temp_classes.py")
        pass
    intermediate_file6 = base_pred + "_temp_classes_table.tsv_temp_DR_class.tsv"
    intermediate_file7 = base_pred + "_temp_classes_table.tsv_temp_IS_class.tsv"
    intermediate_file8 = base_pred + "_temp_classes_table.tsv_temp_classes_mirnas.txt"
    try:
        os.system("python3 make_final_classes.py " + intermediate_file6 + " " + intermediate_file7 + " " + intermediate_file8 + " " + out)
    except:
        print("Error: Failed to execute python3 make_final_classes.py")
        pass

    # Removes temp files
    remove_temp_files()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make miRF classes helper")
    parser.add_argument("rnafold", type=str, help="RNAfold file name")
    parser.add_argument("fasta", type=str, help="FASTA file name")
    parser.add_argument("base_pred", type=str, help="Annotation file name")
    parser.add_argument("out", type=str, help="Output file name - you can choose a name for the output file")

    # Parsing arguments
    args = parser.parse_args()

    # Check if the output file already exists
    if os.path.exists(args.out):
        response = input(f"Warning: The output file '{args.out}' already exists. Do you want to overwrite it? (y/n): ").strip().lower()
        if response != 'y':
            print("Operation cancelled by the user.")
            sys.exit(1)

    # Passing argvs to the runner function
    runner(args.rnafold, args.fasta, args.base_pred, args.out)
