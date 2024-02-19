# Required packages
library(caret)
library(mlbench)
library(dplyr)
library("tidyr")

                          ################################################################### TRAINING ######################################################################### 

# The training dataset consists of 168 negative examples and 163 positive examples (rows)
# For each row, representative features (columns) have been calculated, able to distiguish between positive and negative entities. 
# A normalization has been performed on counts scaled by the hairpin lengths, while the hairpin lenghts have been divided by 100. Frequences, instead remained as they are in input. 

# MODEL -> the algorithm is random forest, tuned with Cross-Validation and the following hyperparameters: mtry=c(5,6,37,50,55,57,119,120,121,122,123)

# Function to perform normalization
performNormalization <- function(data) {
  normalized_data <- data[, 1:36] / data$`hairpin length`
  normalized_data[, 37] <- data[, 37] / 100
  normalized_data <- cbind(normalized_data, data[, 38:71] / data$`hairpin length`)
  normalized_data <- cbind(normalized_data, data[, 72:81])
  normalized_data <- cbind(normalized_data, data[, 82:124] / data$`hairpin length`)
  normalized_data <- cbind(normalized_data, data[, 125:125])
  colnames(normalized_data)[37] <- 'hairpin length'
  colnames(normalized_data)[125] <- 'real miRNA'
  return(normalized_data)
}

# Training dataframe construction 
tab168_new <- read.csv("features_table_for_only_168.txt", row.names = 1, header=TRUE, sep="\t", check.names = F)
tab163_new <- read.csv("features_table_for_only_163.txt", row.names = 1, header=TRUE, sep="\t", check.names = F)
tab168_new$`real miRNA` <- 'no'
tab163_new$`real miRNA` <- 'yes'
tab331_new <- rbind(tab168_new, tab163_new)
column_index331_new <- which(names(tab331_new) == "real miRNA")
new_column_order331_new <- c(names(tab331_new)[-column_index331_new], "real miRNA")
df331_new <- tab331_new[new_column_order331_new]

# Normalization using the function
miR_normalized_df <- performNormalization(df331_new)


# MODEL
set.seed(825)
fitC331_new <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
set.seed(825)
miR_MODEL <- train(`real miRNA` ~ ., data = miR_normalized_df, method = "rf", trControl = fitC331_new, tuneGrid = expand.grid(mtry=c(5,6,37,50,55,57,119,120,121,122,123)))                         

# saveRDS(miR_MODEL, file = "trained_model_new.RDS")
