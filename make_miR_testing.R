                          ################################################################### TESTING ######################################################################### 

# I have decided to test the METHOD on 2 independent datasets:
# TEST 1 consists on 1000 negative esamples and 394 positive examples 
# TEST 2 consistes on 193 negative examples and 193 positive examples

                                                                          #################### TEST 1 #######################

# TEST1 dataframe construction 
hseq_ok = read.csv("features_table_for_hseqCORRECT.out",row.names = 1, header=TRUE, sep="\t", check.names = F)
miRNA_names <- sub(" .*", "", rownames(hseq_ok)) # Extract miRNA names from row names of hseq_tab
hseq_ok$miRNA_name <- miRNA_names
examples_30_plus_163 = read.csv("30_plus_163_mirna.txt",header=F)
mirgene_names = read.csv("mirgene_names_uniq.txt", header=F)
examples_30_plus_163_names <- examples_30_plus_163$V1
mirgene_names <- mirgene_names$V1
unique_mirnas <- mirgene_names[!mirgene_names %in% examples_30_plus_163_names] # Elements in mirgene_names not in examples_30_plus_163_names
unique_mirnas = as.data.frame(unique_mirnas)
names_to_extract_mirgene <- unique_mirnas$unique_mirnas
subset_rows_mirgene <- hseq_ok[hseq_ok$miRNA_name %in% names_to_extract_mirgene, ]
subset_rows_mirgene = subset_rows_mirgene[,1:125]
subset_rows_mirgene$`real miRNA`="yes"
tab30 = read.csv("features_table_for_30_hsa_mir", row.names = 1, header=TRUE, sep="\t", check.names = F)
tab30$`real miRNA`="yes"
mirgene_and_30_mirnas = rbind(tab30,subset_rows_mirgene)
tab1000 = read.csv("features_table_for_1000_cds", row.names = 1, header=TRUE, sep="\t", check.names = F)
tab1000$`real miRNA`="no"
TeS_1 = rbind(tab1000, mirgene_and_30_mirnas) # This is the new test (1394 examples (1000 no and 394 yes))
column_index_test2 <- which(names(TeS_1) == "real miRNA")
new_column_order_test1 <- c(names(TeS_1)[-column_index_test2], "real miRNA")
df_test1 <- TeS_1[new_column_order_test1]

# Normalization
miR_Radial_test1 <- performNormalization(df_test1)

# Testing
set.seed(825)
test1_testingMODEL <- predict(miR_MODEL, miR_Radial_test1)
cm_test1 <- confusionMatrix(positive="yes", test1_testingMODEL, as.factor(miR_Radial_test1$`real miRNA`))

                       
                                                                      #################### TEST 2 #######################

# TEST2 dataframe construction 
tab_pseudo = read.csv("features_table_for_193_pseudo.txt", row.names = 1, header=TRUE, sep="\t", check.names = F)
names_cons = read.csv("mirna_conserved_names.txt",header = F)
miRNA_names <- sub(" .*", "", rownames(hseq_ok)) # Extract miRNA names from row names of hseq_tab
hseq_ok$miRNA_name <- miRNA_names
names_to_extract <- names_cons$V1
subset_rows <- hseq_ok[hseq_ok$miRNA_name %in% names_to_extract, ]
tab_193 = subset_rows[,1:125]
tab_pseudo$`real miRNA`='no'
tab_193$`real miRNA`='yes'
df_pseudo_conserved = rbind(tab_pseudo,tab_193)
column_indexpc <- which(names(df_pseudo_conserved) == "real miRNA")
new_column_orderpc <- c(names(df_pseudo_conserved)[-column_indexpc], "real miRNA")
df_test2 <- df_pseudo_conserved[new_column_orderpc]

# Normalization
miR_Radial_test2 <- performNormalization(df_test2)

# Testing
set.seed(825)
test2_testingMODEL <- predict(miR_MODEL, miR_Radial_test2)
cm_test2 <- confusionMatrix(positive="yes", test2_testingMODEL, as.factor(miR_Radial_test2$`real miRNA`))


                                                                     #################### RESULTS #######################


