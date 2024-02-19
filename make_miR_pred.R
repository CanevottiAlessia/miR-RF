#!/usr/bin/env Rscript

# Check if the trained model exists, if not, train and save it
if (!file.exists("trained_model.RDS")) {
  train_model_and_save()
}

args=commandArgs(trailingOnly=TRUE)
name_table=args[1]
name_final=args[2]
trained_model <- readRDS("trained_model_new.RDS")

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
input_R[is.na(input_R)] <- 0

column_index <- which(names(input_R) == "real miRNA")
new_column_order <- c(names(input_R)[-column_index], "real miRNA")
df_1917 <-input_R[new_column_order]
new_1917_norm_by_len = df_1917[,1:124]/df_1917$`hairpin length`
new_1917_normL <- cbind(new_1917_norm_by_len, df_1917[, 125:125])
colnames(new_1917_normL)[125] <- 'real miRNA'
set.seed(825)

# Use the loaded trained model to make predictions
output <- predict_output(trained_model, new_1917_normL)
colnames(output) <- c("miRNA name", "prediction")

# Write the updated output to a file with the specified column names
write.table(output, name_final, sep = "\t", col.names = TRUE, row.names = FALSE)

#write.table(output, "final_output.tsv")
