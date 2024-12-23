import sys
import os

def runner(inp):
    try:
        os.system("python3 PY_miR_features_extraction.py " + inp)
    except:
        print("Error: Failed to execute. Please, enter the input file name.")
        pass
    intermediate_file="features_table_for_" + inp
    try:
        os.system("Rscript --vanilla make_miR_pred_classes2.R " + intermediate_file) # + " " + out)
    except:
        print("Error: Failed to execute. Please, enter the output file name.")
        pass
    try:
        os.remove(intermediate_file)
    except OSError as e:
        print(f"Error removing the file: {e}")



if __name__=="__main__":
     user_file = sys.argv[1]
     runner(user_file)
