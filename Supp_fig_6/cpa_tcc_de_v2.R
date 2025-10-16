#read metadata /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/TCC/tcc_h5s_0522/final_df.csv
meta_data = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/tcc_h5s_0522/final_df.csv")
#in meta_data only keep rows where cell_type is not Blood, Endothelial, Connective tissue or Immune

print(unique(meta_data$cell_type))
filenames = meta_data$filename


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("devtools")    # only if devtools not yet installed
#BiocManager::install("biomaRt")
#BiocManager::install("pachterlab/sleuth",force=TRUE)

library('sleuth')
library(biomaRt)

filenames = paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/tcc_h5s_0522/", meta_data$filename)
#attempt to load in all h5 files, remove ones that cannot be loaded



metadata = data.frame(cell_type = meta_data$cell_type, 
                      path = filenames,
                      sample = paste0(meta_data$SRA_ID, meta_data$cell_type))





#rename columns to tx2gene to target_id and ens_gene
library(biomaRt)
mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")


ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)


ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
head(ttg)

#remove cols in metadata that are not entirely NAN free
#interate thru cols
for (col in colnames(metadata)){
  if (sum(is.na(metadata[,col])) > 0){
    metadata[,col] <- NULL
  }
}



cell_types = unique(metadata$cell_type)


library(tidyverse)
new_metadata = metadata
for (ct in cell_types){
  new_metadata = new_metadata %>% mutate(ct = ifelse(cell_type == ct, ct, "rest"))
  #rename this column to ct value
  colnames(new_metadata)[ncol(new_metadata)] = ct
}
  




all_results_isoform = list()

grouping_1 <- list(list("Stem"), list("Melanotrophs","Corticotrophs",
                                            "Somatotrophs", "Lactotrophs","Thyrotrophs",
                                            "Gonadotrophs"
))

grouping_2 <- list(list("Gonadotrophs"), list("Melanotrophs","Corticotrophs",
                                              "Somatotrophs", "Lactotrophs","Thyrotrophs",
                                              "Stem"
                                              
))

grouping_3 <- list(list("Melanotrophs","Corticotrophs"), list(
  "Somatotrophs", "Lactotrophs","Thyrotrophs",
  "Gonadotrophs", "Stem"
))

grouping_4 <- list(list("Melanotrophs"), list("Corticotrophs"))

grouping_5 <- list(list("Somatotrophs", "Lactotrophs","Thyrotrophs"), list(
  "Gonadotrophs", "Stem", "Melanotrophs",
  "Corticotrophs"))

grouping_6 <- list(list("Lactotrophs"), list("Somatotrophs", "Thyrotrophs"))

grouping_7 <- list(list("Somatotrophs"), list("Lactotrophs", "Thyrotrophs"))

grouping_8 <- list(list("Thyrotrophs"), list("Lactotrophs", "Somatotrophs"))

# List of all groupings
all_groupings <- list(grouping_1, grouping_2, grouping_3, grouping_4, 
                      grouping_5, grouping_6, grouping_7, grouping_8)






find_significant_genes <- function(metadata, grouping1, grouping2) {
  #define new_metadata, where make a grouping variable, and assign one if cell type in grouping1
  print("***************")
  print("Running")
  print("***************")
  new_metadata = metadata
  #only keep those cell types that are either in grouping1 or grouping2
  new_metadata = new_metadata[new_metadata$cell_type %in% c(grouping1, grouping2),]
  
  new_metadata = new_metadata %>% mutate(grouping = ifelse(cell_type %in% grouping1, 1, 0))
  
  so_orig <- sleuth_prep(new_metadata, target_mapping = ttg,max_bootstrap=10,
                         aggregation_column = 'ens_gene', extra_bootstrap_summary = FALSE,
                         transform_fun_counts = function(x) log2(x + 0.5), gene_mode = FALSE)
  
  
  so <- sleuth_fit(so_orig, ~grouping, 'full')
  
  so <- sleuth_wt(so, "grouping")
  sleuth_table_gene <- sleuth_results(so, "grouping",pval_aggregate = FALSE)
  return(sleuth_table_gene)
}
  
  
process_grouping <- function(metadata, grouping, grouping_name) {
  result <- find_significant_genes(metadata, grouping[[1]], grouping[[2]])
  if (nrow(result) > 0) {
    result$grouping <- grouping_name
    return(result)
  }
}

# Process all groupings
results_list <- lapply(seq_along(all_groupings), function(i) {
  process_grouping(metadata,all_groupings[[i]], paste0("grouping_", i))
})


#save results_list to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/TCC/results_lineage_isoform.rds
saveRDS(results_list, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/results_lineage_isoform_0728.rds")
results_list <- readRDS("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/results_lineage_isoform_0728.rds")


res1 <- results_list[[1]]
#save as csv
write.csv(res1, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/results_lineage_isoform_grouping_1.csv")

find_significant_genes <- function(metadata, grouping1, grouping2) {
  #define new_metadata, where make a grouping variable, and assign one if cell type in grouping1
  print("***************")
  print("Running")
  print("***************")
  new_metadata = metadata
  #only keep those cell types that are either in grouping1 or grouping2
  new_metadata = new_metadata[new_metadata$cell_type %in% c(grouping1, grouping2),]
  
  new_metadata = new_metadata %>% mutate(grouping = ifelse(cell_type %in% grouping1, 1, 0))
  
  so_orig <- sleuth_prep(new_metadata, target_mapping = ttg,max_bootstrap=10,
                         aggregation_column = 'ens_gene', extra_bootstrap_summary = FALSE,
                         transform_fun_counts = function(x) log2(x + 0.5), gene_mode = TRUE)
  
  
  so <- sleuth_fit(so_orig, ~grouping, 'full')
  
  so <- sleuth_wt(so, "grouping")
  sleuth_table_gene <- sleuth_results(so, "grouping",pval_aggregate = FALSE)
  return(sleuth_table_gene)
}


process_grouping <- function(metadata, grouping, grouping_name) {
  result <- find_significant_genes(metadata, grouping[[1]], grouping[[2]])
  if (nrow(result) > 0) {
    result$grouping <- grouping_name
    return(result)
  }
}

# Process all groupings
results_list_gene <- lapply(seq_along(all_groupings), function(i) {
  process_grouping(metadata,all_groupings[[i]], paste0("grouping_", i))
})

results_list_gene[[1]]


#save 
saveRDS(results_list_gene, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/results_lineage_gene_0728.rds")


library(tidyverse)


concatted_for_each= list()
for (i in 1:length(results_list)){
  #modify the gene_mode df by addig an isoform_qval col
  isoform_df =  results_list[[i]]
  gene_df = results_list_gene[[i]]
  gene_qvals = c()
  #iterate through ens_gene
  for (ens_gene in isoform_df$ens_gene){
    #get all isoform qvals for that gene
    gene_qval = gene_df[gene_df$target_id == ens_gene,]$qval
    if (is.na(ens_gene)){
      gene_qval = NA
    }
    if (length(gene_qval) == 0){
      gene_qval = NA }
  
    gene_qvals = c(gene_qvals, gene_qval[1])
  }
  isoform_df$gene_qval = gene_qvals
  #add isoform_df$n_target_id which is number of times given target_id appears in the isoform_df
  isoform_df <- as.data.frame(isoform_df)
  
  isoform_df = isoform_df %>% group_by(ext_gene) %>% mutate(n_target_id = n())
  concatted_for_each[[i]] = isoform_df
}

tab <-concatted_for_each[[1]]


#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/TCC/
saveRDS(concatted_for_each, file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/concatted_for_each_0728.rds")


#filter for cases where qval < 0.05 and 
finalised_results = list()
for (i in 1:length(concatted_for_each)){
  df = concatted_for_each[[i]]
  df = df[df$qval < 0.05,]
  #and also isoform qval smaller than gene qval
  df = df[df$qval < df$gene_qval,]
  finalised_results[[i]] = df
}


#bind these into a single one and add col saying which coef
finalised_results_final = list()
for (i in 1:length(finalised_results)){
  df = finalised_results[[i]]
  df$grouping = c("Grouping_1", "Grouping_2", "Grouping_3", "Grouping_4", "Grouping_5", "Grouping_6", "Grouping_7", "Grouping_8")[i]
  finalised_results_final[[i]] = df
}
#bind
finalised_results_final = do.call(rbind, finalised_results_final)
#add column called log10_diff_in_qvals
finalised_results_final$log10_diff_in_qvals = log10(finalised_results_final$qval) - log10(finalised_results_final$gene_qval)
#order
#save

saveRDS(finalised_results_final , file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/finalised_results_final_0728.rds")
finalised_results_final <- readRDS("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/finalised_results_final_0728.rds")


finalised_results_final
#remove rows with nans
finalised_results_final = finalised_results_final[!is.na(finalised_results_final$ext_gene),]
#remove where n_target_id is nan
finalised_results_final = finalised_results_final[!is.na(finalised_results_final$n_target_id),]
#for each ens_gene, only keep the row with the lowest log10_diff_in_qvals
finalised_results_final = finalised_results_final[order(finalised_results_final$log10_diff_in_qvals),]
finalised_results_final = finalised_results_final[!duplicated(finalised_results_final$ens_gene),]
#keep those where log10_diff_in_qvals is less than -1
finalised_results_final = finalised_results_final[finalised_results_final$log10_diff_in_qvals < -1,]

#read in those isoforms that we have uniquely mapping reads onto 
#/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/isoforms/v_0.01/filtered_transcripts_list.csv
isoforms = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/isoforms/v_0.01/filtered_transcripts_list.csv")
#keep only finalised_results_final where target_id is in isoforms
#reset index of finalised_results_final
rownames(finalised_results_final) = 1:nrow(finalised_results_final)
finalised_results_final = finalised_results_final[finalised_results_final$target_id %in% isoforms$X0,]
finalised_results_final

#save as csv called isoforms results to plot
write.csv(finalised_results_final, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/isoforms_results_to_plot.csv")

finalised_results_final = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/isoforms_results_to_plot.csv")


#keep only tfs
tf_list <- read.table("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/allTFs_mm.txt", header = FALSE, stringsAsFactors = FALSE)
tf_list <- tf_list$V1

#keep where ext_gene is in tf_list
finalised_results_final = finalised_results_final[finalised_results_final$ext_gene %in% tf_list,]

finalised_results_final












###################
#finally, a comparison between the pituitary lineage and the rest of the cell types
###################




all_results_isoform = list()

grouping_9 <- list(list("Stem","Melanotrophs","Corticotrophs",
                        "Somatotrophs", "Lactotrophs","Thyrotrophs",
                        "Gonadotrophs"), list("Endothelial","Pituicytes","Immune",
                                              "Mesenchymal"))

# List of all groupings
all_groupings <- list(grouping_9)






find_significant_genes <- function(metadata, grouping1, grouping2) {
  #define new_metadata, where make a grouping variable, and assign one if cell type in grouping1
  print("***************")
  print("Running")
  print("***************")
  new_metadata = metadata
  #only keep those cell types that are either in grouping1 or grouping2
  new_metadata = new_metadata[new_metadata$cell_type %in% c(grouping1, grouping2),]
  
  new_metadata = new_metadata %>% mutate(grouping = ifelse(cell_type %in% grouping1, 1, 0))
  
  so_orig <- sleuth_prep(new_metadata, target_mapping = ttg,max_bootstrap=10,
                         aggregation_column = 'ens_gene', extra_bootstrap_summary = FALSE,
                         transform_fun_counts = function(x) log2(x + 0.5), gene_mode = FALSE)
  
  
  so <- sleuth_fit(so_orig, ~grouping, 'full')
  
  so <- sleuth_wt(so, "grouping")
  sleuth_table_gene <- sleuth_results(so, "grouping",pval_aggregate = FALSE)
  return(sleuth_table_gene)
}


process_grouping <- function(metadata, grouping, grouping_name) {
  result <- find_significant_genes(metadata, grouping[[1]], grouping[[2]])
  if (nrow(result) > 0) {
    result$grouping <- grouping_name
    return(result)
  }
}

# Process all groupings
results_list <- lapply(seq_along(all_groupings), function(i) {
  process_grouping(metadata,all_groupings[[i]], paste0("grouping_", i))
})


#save results_list to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/TCC/results_lineage_isoform.rds
saveRDS(results_list, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/results_pit_vs_nopit_isoform_0728.rds")
results_list <- readRDS("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/results_pit_vs_nopit_isoform_0728.rds")








find_significant_genes <- function(metadata, grouping1, grouping2) {
  #define new_metadata, where make a grouping variable, and assign one if cell type in grouping1
  print("***************")
  print("Running")
  print("***************")
  new_metadata = metadata
  #only keep those cell types that are either in grouping1 or grouping2
  new_metadata = new_metadata[new_metadata$cell_type %in% c(grouping1, grouping2),]
  
  new_metadata = new_metadata %>% mutate(grouping = ifelse(cell_type %in% grouping1, 1, 0))
  
  so_orig <- sleuth_prep(new_metadata, target_mapping = ttg,max_bootstrap=10,
                         aggregation_column = 'ens_gene', extra_bootstrap_summary = FALSE,
                         transform_fun_counts = function(x) log2(x + 0.5), gene_mode = TRUE)
  
  
  so <- sleuth_fit(so_orig, ~grouping, 'full')
  
  so <- sleuth_wt(so, "grouping")
  sleuth_table_gene <- sleuth_results(so, "grouping",pval_aggregate = FALSE)
  return(sleuth_table_gene)
}


process_grouping <- function(metadata, grouping, grouping_name) {
  result <- find_significant_genes(metadata, grouping[[1]], grouping[[2]])
  if (nrow(result) > 0) {
    result$grouping <- grouping_name
    return(result)
  }
}

# Process all groupings
results_list_gene <- lapply(seq_along(all_groupings), function(i) {
  process_grouping(metadata,all_groupings[[i]], paste0("grouping_", i))
})

results_list_gene[[1]]


#save 
saveRDS(results_list_gene, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/results_pit_vs_nopit_gene_0728.rds")


library(tidyverse)


concatted_for_each= list()
for (i in 1:length(results_list)){
  #modify the gene_mode df by addig an isoform_qval col
  isoform_df =  results_list[[i]]
  gene_df = results_list_gene[[i]]
  gene_qvals = c()
  #iterate through ens_gene
  for (ens_gene in isoform_df$ens_gene){
    #get all isoform qvals for that gene
    gene_qval = gene_df[gene_df$target_id == ens_gene,]$qval
    if (is.na(ens_gene)){
      gene_qval = NA
    }
    if (length(gene_qval) == 0){
      gene_qval = NA }
    
    gene_qvals = c(gene_qvals, gene_qval[1])
  }
  isoform_df$gene_qval = gene_qvals
  #add isoform_df$n_target_id which is number of times given target_id appears in the isoform_df
  isoform_df <- as.data.frame(isoform_df)
  
  isoform_df = isoform_df %>% group_by(ext_gene) %>% mutate(n_target_id = n())
  concatted_for_each[[i]] = isoform_df
}

tab <-concatted_for_each[[1]]


#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/TCC/
saveRDS(concatted_for_each, file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/concatted_for_pit_vs_nopit_0728.rds")


#filter for cases where qval < 0.05 and 
finalised_results = list()
for (i in 1:length(concatted_for_each)){
  df = concatted_for_each[[i]]
  df = df[df$qval < 0.05,]
  #and also isoform qval smaller than gene qval
  df = df[df$qval < df$gene_qval,]
  finalised_results[[i]] = df
}


#bind these into a single one and add col saying which coef
finalised_results_final = list()
for (i in 1:length(finalised_results)){
  df = finalised_results[[i]]
  df$grouping = c("Grouping_9")[i]
  finalised_results_final[[i]] = df
}
#bind
finalised_results_final = do.call(rbind, finalised_results_final)
#add column called log10_diff_in_qvals
finalised_results_final$log10_diff_in_qvals = log10(finalised_results_final$qval) - log10(finalised_results_final$gene_qval)
#order
#save

saveRDS(finalised_results_final , file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tcc/finalised_pit_vs_nopit_results_final_0728.rds")







































































write.table(sleuth_table_gene, file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_one/Proposal/pituitary atlas/kallisto_results_gene_aggr_russell.tsv", sep = "\t", quote = FALSE)

sleuth_table_transcript <- sleuth_results(so, 'conditionb',pval_aggregat=FALSE)
write.table(sleuth_table_transcript, file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_one/Proposal/pituitary atlas/kallisto_results_tx_russell.tsv", sep = "\t", quote = FALSE)


so <- sleuth_prep(metadata, target_mapping = ttg,max_bootstrap=100,
                  aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE,
                  transform_fun_counts = function(x) log2(x + 0.5),gene_mode = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_wt(so, 'conditionb')
sleuth_table_gene_gene <- sleuth_results(so, 'conditionb')

write.table(sleuth_table_gene_gene, file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_one/Proposal/pituitary atlas/kallisto_results_genemode_russell.tsv", sep = "\t", quote = FALSE)


