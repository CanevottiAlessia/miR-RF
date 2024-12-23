import pandas as pd
import sys
import os

if __name__=="__main__":
   user_file = sys.argv[1]
   user_file2 = sys.argv[2]

# BASE

file_path = user_file 
df1 = pd.read_csv(file_path, sep='\t')
df1 = df1[['miRNA name', 'prediction']]

# POTENTIAL

file_path = user_file2 
df = pd.read_csv(file_path, sep='\t')

# Crea nuove colonne per mutation e position
df['mutation'] = df['miRNA name'].str.extract(r'([A-Z]-[A-Z])')  # Estracts U-A, U-C, ecc.
df['position'] = df['miRNA name'].str.extract(r'(pos=\d+)')       # Estracts pos=numero

# Rimuove mutation e position dalla colonna 'miRNA name'
df['miRNA name'] = df['miRNA name'].str.replace(r' [A-Z]-[A-Z], pos=\d+', '', regex=True)

# Mantieni solo le colonne richieste
df = df[['miRNA name', 'mutation', 'position', 'prediction']]

# MERGE THE TABLES

df1['prediction'] = df1['prediction'].replace(2, "YES")
df1['prediction'] = df1['prediction'].replace(1, "NO")

df['prediction'] = df['prediction'].replace(2, "YES")
df['prediction'] = df['prediction'].replace(1, "NO")

merging_df = pd.merge(df1, df, on='miRNA name')
#merging_df.to_csv("merged_tables.tsv", sep="\t", index=False)

# DEFINE IF THE MUTATION IS DISRUPTING OR NOT (miRNA status)

def get_mirna_status(row):
    if row["prediction_x"] == "YES" and row["prediction_y"] == "YES":
        return "neutral"
    elif row["prediction_x"] == "YES" and row["prediction_y"] == "NO":
        return "deactivated"
    elif row["prediction_x"] == "NO" and row["prediction_y"] == "NO":
        return "noImpact"
    elif row["prediction_x"] == "NO" and row["prediction_y"] == "YES":
        return "activated"
merging_df["Mirna_Status"] = merging_df.apply(get_mirna_status, axis=1)

merging_df.to_csv(user_file + "_temp_classes_table.tsv", sep="\t", index=False) #merging_df.to_csv("temp_classes_table.tsv", sep="\t", index=False)
