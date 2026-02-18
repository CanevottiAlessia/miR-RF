#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
user_file  <- args[1]
out_folder <- args[2]

# se arg1 è già la features table, usala così com'è
if (grepl("_features_table\\.tsv$", basename(user_file))) {
  name_table <- user_file
  base_name <- sub("_features_table\\.tsv$", "", basename(user_file))
} else {
  # altrimenti, costruiscila come faceva il Python
  name_table <- file.path(out_folder, paste0(basename(user_file), "_features_table.tsv"))
  base_name <- basename(user_file)
}

# output sempre dentro out_folder
name_final <- file.path(out_folder, paste0(base_name, "_pred.tsv"))

# 1) Check che la feature table esista
if (!file.exists(name_table)) {
  stop(paste("Input feature table not found:", name_table))
}

# 2) Check che il modello esista nella stessa cartella dello script R
model_file <- "trained_model_new.RDS"

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
script_path <- sub("^--file=", "", file_arg[1])
script_dir <- dirname(normalizePath(script_path))

model_path <- file.path(script_dir, model_file)

if (!file.exists(model_path)) {
  stop(paste("Trained model not found:", model_path, "\n"))
}
trained_model <- readRDS(model_path)

# PREDICTION FUNCTION
predict_output <- function(trained_model, new_data) {
  validating_rf <- predict(trained_model, new_data)
  output <- cbind(rownames(new_data), validating_rf)
  return(output)
}

input_R <- read.table(name_table, row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
for (col in colnames(input_R)) {
  input_R[, col] <- as.numeric(input_R[, col])
}

# .................................... MANAGEMENT ..........................................

# Management of too long hairpins (>5 loops)
pseudo_count <- 1e-6
input_R[is.na(input_R)] <- pseudo_count
# Management of null hairpins (miRNAs without loops after the SNP substitution)
input_R[input_R$`hairpin length` == 0, ] <- lapply(input_R[input_R$`hairpin length` == 0, ],
                                                           function(x) ifelse(x == 0, 1e-6, x))
# I need to temporaly assign 0 to "real miRNA" column
input_R$`real miRNA`=0

# .................................... PROCEED ............................................
column_index <- which(names(input_R) == "real miRNA")
new_column_order <- c(names(input_R)[-column_index], "real miRNA")
df_1917 <- input_R[new_column_order]

# Normalization
performNormalization <- function(data) {
  normalized_data <- data[, 1:36] / data$`hairpin length`
  normalized_data[, 37] <- data[, 37] / 100
  normalized_data <- cbind(normalized_data, data[, 38:71] / data$`hairpin length`)
  normalized_data <- cbind(normalized_data, data[, 72:81])
  normalized_data <- cbind(normalized_data, data[, 82:124] / data$`hairpin length`)
  normalized_data <- cbind(normalized_data, data[, 125:125])
  colnames(normalized_data)[37] <- "hairpin length"
  colnames(normalized_data)[125] <- "real miRNA"
  return(normalized_data)
}

new_1917_normL <- performNormalization(df_1917)

# Predict
set.seed(825)
output <- predict_output(trained_model, new_1917_normL)
colnames(output) <- c("miRNA name", "prediction")

# Write output
write.table(output, name_final, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
