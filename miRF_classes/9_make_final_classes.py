import pandas as pd
import sys
import os

if __name__=="__main__":
    user_file = sys.argv[1]
    user_file2 = sys.argv[2]
    user_file3 = sys.argv[3]
    user_file4 = sys.argv[4]

# Read the input files
is_class = pd.read_csv(user_file2, sep="\t") #"IS_class.tsv", sep="\t")
is_class = is_class.iloc[:, [0, 6]]  # Take the first and the last columns
is_class.columns = ["mirna_name", "status"]

dr_class = pd.read_csv(user_file, sep="\t") #"DR_class.tsv", sep="\t")
dr_class = dr_class.iloc[:, [0, 6]]  # Take the first and the last columns
dr_class.columns = ["mirna_name", "status"]

classes_temp = pd.read_csv(user_file3, sep="\t", header=None, names=['mirna_name', 'status']) 
# Rename classes
classes_temp["status"] = classes_temp["status"].replace({
    "Dispensable": "D",
    "Inducible": "I",
    "Spurious": "S",
    "Resilient": "R"
})

# Combine IS and DR classes
tab_ISDR_temp = pd.concat([is_class, dr_class], ignore_index=True)
tab_ISDR_temp = tab_ISDR_temp[tab_ISDR_temp["status"].isin(["IS", "DR"])]

# Join and update statuses
class_temp_updated = pd.merge(classes_temp, tab_ISDR_temp, on="mirna_name", how="left", suffixes=("_classes_temp", "_tab_ISDR"))
class_temp_updated["status"] = class_temp_updated["status_tab_ISDR"].combine_first(class_temp_updated["status_classes_temp"])

# Reclassify statuses
class_temp_updated["new_status"] = class_temp_updated["status"].map({
    "D": "D",
    "R": "R",
    "DR": "R",
    "S": "S",
    "IS": "S",
    "I": "I"
})

# Select relevant columns
new_classes_final = class_temp_updated[["mirna_name", "new_status"]]

# Write the output to a file
new_classes_final.to_csv(user_file4 + "_final_classes.tsv", index=False, sep="\t")
