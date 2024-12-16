import sys
import os

def runner(rnafold,fasta,base_pred):
    try:
        os.system("python3 3_get_fasta_for_insertion.py " + rnafold)
    except:
        print("Error: Failed to execute. Please, enter the input file name.")
        pass
    intermediate_file = rnafold + "_final_hairpins.txt"
    try:
        os.system("python3 4_SNP_insertion.py " + fasta + " " + intermediate_file)
    except:
        print("Error: Failed to execute python3 4_SNP_insertion.py")
        pass
    intermediate_file2 = intermediate_file + "_potential_RNAfold.txt"
    try:
        os.system("python3 2_miR_application.py " + intermediate_file2)
    except:
        print("Error: Failed to execute python3 2_miR_application.py")
        pass
    intermediate_file3 = "features_table_for_" + intermediate_file2 + "_temp_pred_potential.tsv"
    try:
        os.system("python3 5_merge_table.py " + base_pred + " " + intermediate_file3)
    except:
        print("Error: Failed to execute python3 5_merge_table.py")
        pass
    intermediate_file4 = base_pred + "_temp_classes_table.tsv"
    try:
        os.system("python3 6_get_lens.py " + rnafold)
    except OSError as e:
        print("Error: Failed to execute python3 6_get_lens.py")
        pass
    intermediate_file5 = rnafold + "_temp_len_hairpins.txt"
    try:
        os.system("python3 7_make_fisher_test.py " + intermediate_file4 + " " + intermediate_file5)
    except:
        print("Error: Failed to execute python3 7_make_fisher_test.py")
        pass
    try:
        os.system("python3 8_make_temp_classes.py " + intermediate_file4)
    except:
        print("Error: Failed to execute python3 8_make_temp_classes.py")
        pass
    intermediate_file6 = base_pred + "_temp_classes_table.tsv_temp_DR_class.tsv"
    intermediate_file7 = base_pred + "_temp_classes_table.tsv_temp_IS_class.tsv"
    intermediate_file8 = base_pred + "_temp_classes_table.tsv_temp_classes_mirnas.txt"
    try:
        os.system("python3 9_make_final_classes.py " + intermediate_file6 + " " + intermediate_file7 + " " + intermediate_file8 + " " + rnafold)
    except:
        print("Error: Failed to execute python3 9_make_final_classes.py")
        pass


if __name__ == "__main__":
    rnafold = sys.argv[1]
    fasta = sys.argv[2]
    base_pred = sys.argv[3]
    runner(rnafold, fasta, base_pred)
