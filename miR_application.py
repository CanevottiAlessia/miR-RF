import sys
import os

def runner(inp, out):
    try:
        os.system("python3 PY_miR_features_extraction.py " + inp)
    except:
        print("Error: Failed to execute. Please, enter the input file name.")
        pass
    intermediate_file="features_table_for_" + inp
    try:
        os.system("Rscript make_pred.R " + intermediate_file + " " + out) 
    except:
        print("Error: Failed to execute. Please, enter the output file name.")
        pass
    try:
        os.remove(intermediate_file)
    except OSError as e:
        print(f"Error removing the file: {e}")



if __name__=="__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 scripts.py <INPUT> <OUTPUT>")
    else:
        user_file = sys.argv[1]
        final_file = sys.argv[2]
        runner(user_file, final_file)
