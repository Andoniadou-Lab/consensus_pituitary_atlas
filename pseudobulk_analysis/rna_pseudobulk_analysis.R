library(DESeq2)
library(Seurat)
library(tidyverse)
library(reticulate)
version="v_0.01"
# Specify the path to the .h5ad file
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/pdatas0828.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

# Convert the data to Seurat object
# Seurat expects the count matrix to be in column-major order (genes x cells)
counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = counts, meta.data = meta_data)
#print unique Author
rownames(seurat_object) <- vars
seurat_object$assignments <- seurat_object$new_cell_type

#normalize
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1000000)
                              
root = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/"
figs_folder = paste(root,"Figures/",sep="")


meta_data = seurat_object@meta.data
print(length(unique(meta_data$Author)))
print(length(unique(meta_data$SRA_ID)))

#print sum psbulk_n_cells
print(sum(meta_data$psbulk_n_cells))
expression_table = seurat_object@assays$RNA$counts
unique(meta_data$SRA_ID)


###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$Modality)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)


#remove where assignment is Blood
expression_table <- expression_table[,meta_data$assignments != "Erythrocytes"]
meta_data <- meta_data[meta_data$assignments != "Erythrocytes",]

#Where Protocol contains Parse set exptype to parse
meta_data$exptype[meta_data$`10X version` == "PARSE_WT_MegaV2"] <- "parse"
meta_data$exptype[meta_data$`10X version` == "Parse_WT"] <- "parse"

unique(meta_data$exptype)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$exptype)
meta_data$Sex <- make.names(meta_data$Comp)

# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)
meta_data$Sex <- gsub("[^A-Za-z0-9_]", "", meta_data$Sex)

#order meta_data based on exptype
#reorder expression_table based on meta_data index
expression_table <- expression_table[,order(meta_data$exptype)]
meta_data <- meta_data[order(meta_data$exptype),]

#remove where age_numeric < 0
expression_table <- expression_table[,meta_data$Age_numeric >= 0]
meta_data <- meta_data[meta_data$Age_numeric >= 0,]
#print all values
print(meta_data$exptype)
#add a new column called index
meta_data$index <- 1:nrow(meta_data)
#make index the rownames
rownames(meta_data) <- meta_data$index


print(length(unique(meta_data$SRA_ID)))


meta_data %>%
  #filter where exptype is not parse
  dplyr::filter(exptype != "parse") %>%
  dplyr::filter(Age_numeric > 10) %>%
  dplyr::filter(Age_numeric < 150) %>%
  dplyr::distinct(SRA_ID, Comp_sex) %>%  
  #filter where Age_numeric >0
  dplyr::count(Comp_sex, name = "num_SRA_IDs")



meta_data %>%
  #filter where exptype is not parse
  dplyr::filter(exptype == "sc") %>%
  dplyr::filter(Age_numeric > 0) %>%
  dplyr::distinct(SRA_ID,exptype) %>%  
  #filter where Age_numeric >0
  dplyr::count(exptype, name = "num_SRA_IDs")


meta_data %>%
  #filter where exptype is not parse
  dplyr::filter(Age_numeric > 0) %>%
  dplyr::distinct(SRA_ID,exptype) %>%  
  #filter where Age_numeric >0
  dplyr::count(exptype, name = "num_SRA_IDs")


expression_table <- expression_table / meta_data$psbulk_n_cells
expression_table <- expression_table * median(meta_data$psbulk_n_cells)

exp_table_sc = expression_table[,meta_data$exptype == "sc"]
meta_data_sc = meta_data[meta_data$exptype == "sc",]

exp_table_sn = expression_table[,meta_data$exptype == "sn"]
meta_data_sn = meta_data[meta_data$exptype == "sn",]

exp_table_multi = expression_table[,meta_data$exptype == "multi_rna"]
meta_data_multi = meta_data[meta_data$exptype == "multi_rna",]

exp_table_parse = expression_table[,meta_data$exptype == "parse"]
meta_data_parse = meta_data[meta_data$exptype == "parse",]

dge = DGEList(counts = expression_table, group = meta_data$assignments)
# Prepare the DGEList object for limma voom
dge_sc <- DGEList(counts = exp_table_sc, group = meta_data_sc$assignments)
dge_sn <- DGEList(counts = exp_table_sn, group = meta_data_sn$assignments)
dge_multi <- DGEList(counts = exp_table_multi, group = meta_data_multi$assignments)
dge_parse <- DGEList(counts = exp_table_parse, group = meta_data_parse$assignments)

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
design <- model.matrix(~ 0 + assignments, data = meta_data)

keep_genes <- filterByExpr(dge, design, min.total.count = 5000, min.prop= 0.5)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)
keep_genes["Prop1"]
dge_sc <- dge_sc[gene_names,]
dge_sn <- dge_sn[gene_names,]
dge_multi <- dge_multi[gene_names,]
dge_parse <- dge_parse[gene_names,]

#print genes names with top counts
rownames(dge)[order(rowSums(dge$counts), decreasing = TRUE)[1:15]]
#also print their total count
rowSums(dge$counts)[order(rowSums(dge$counts), decreasing = TRUE)[1:15]]

lib_sizes_sc = colSums(dge_sc$counts[-order(rowSums(dge_sc$counts), decreasing = TRUE)[1:15],])
lib_sizes_sn = colSums(dge_sn$counts[-order(rowSums(dge_sn$counts), decreasing = TRUE)[1:15],])
lib_sizes_multi = colSums(dge_multi$counts[-order(rowSums(dge_multi$counts), decreasing = TRUE)[1:15],])
lib_sizes_parse = colSums(dge_parse$counts[-order(rowSums(dge_parse$counts), decreasing = TRUE)[1:15],])

#add these lib_size to dge
dge_sc$samples$lib.size = lib_sizes_sc
dge_sn$samples$lib.size = lib_sizes_sn
dge_multi$samples$lib.size = lib_sizes_multi
dge_parse$samples$lib.size = lib_sizes_parse

#norm
dge_sc <- calcNormFactors(dge_sc, method = "TMM")
dge_sn <- calcNormFactors(dge_sn, method = "TMM")
dge_multi <- calcNormFactors(dge_multi, method = "TMM")
dge_parse <- calcNormFactors(dge_parse, method = "TMM")

#uniting these library sizes and norm factors in dge
dge$samples = rbind(dge_multi$samples,dge_parse$samples, dge_sc$samples,dge_sn$samples)

meta_data = rbind(meta_data_multi, meta_data_parse, meta_data_sc,meta_data_sn)

dge_original <- dge
meta_data_original <- meta_data
expression_table_original <- expression_table

# Create the fixed effects design matrix
meta_data$assignments_exptype <- paste(meta_data$assignments, meta_data$exptype, sep = "_")

# Create the design matrix
design <- model.matrix(~ 0 + assignments, data = meta_data) #+ Sex + assignments:Age_numeric  + assignments:Sex, data = meta_data)
design_sc <- model.matrix(~ 0 + assignments, data = meta_data_sc)
design_sn <- model.matrix(~ 0 + assignments, data = meta_data_sn)
design_multi <- model.matrix(~ 0 + assignments, data = meta_data_multi)
design_parse <- model.matrix(~ 0 + assignments, data = meta_data_parse)

#running voom
fit_sc <- voom(dge_sc, design_sc, plot=TRUE)
fit_sn <- voom(dge_sn, design_sn, plot=TRUE)
fit_multi <- voom(dge_multi, design_multi, plot=TRUE)
fit_parse <- voom(dge_parse, design_parse, plot=TRUE)
fit <- voom(dge, design, plot=TRUE)

#is fit the same as cbind(fit_multi,fit_sc,fit_sn)
fit2=cbind(fit_multi,fit_parse,fit_sc,fit_sn)

#double running this based on
#https://support.bioconductor.org/p/114663/
#first changing the weights to the separate ones
fit$weights <- fit2$weights
corfit <- duplicateCorrelation(fit, design, block=meta_data$exptype)
corfit$consensus.correlation
fit <- voom(dge, design, plot=TRUE, block=meta_data$exptype, correlation=corfit$consensus.correlation)
corfit <- duplicateCorrelation(fit, design, block=meta_data$exptype)
corfit$consensus.correlation

fit <- lmFit(fit, design, block=meta_data$exptype, correlation=corfit$consensus.correlation)


###
# Make column names valid
colnames(design) <- make.names(colnames(design), unique = TRUE)

print(head(coef(fit)))

#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
write.csv(as.data.frame(coef(fit)), "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/coef.csv")
coefs_table <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/coef.csv")
#save to epitome

write.csv(as.data.frame(coef(fit)), "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/heatmap/v_0.01/coef.csv")

#turn coef(fit) into a dataframe
coef_df = as.data.frame(coef(fit))

coef_df$gene = rownames(coef_df)

# Create contrast
contrast <- makeContrasts(
  #assignmentsCorticotrophs - assignmentsMelanotrophs, 
  assignmentsSomatotrophs-assignmentsStem_cells,
  levels = design
)


# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, trend = TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf, sort.by = "P")
top_genes$genes <- rownames(top_genes)
print(head(top_genes))


# Initialize a list to store upregulated genes for each cell type
all_upregulated_genes <- list()


# Get unique cell types
celltypes <- unique(meta_data$assignments)


# Function to compare one cell type against all others
compare_celltype <- function(celltype, other_celltypes, design, fit, log2fc=-Inf) {
  significant_genes <- data.frame()
  
  for (other in other_celltypes) {
    contrast_name <- paste0("assignments", celltype)
    other_contrast_name <- paste0("assignments", other)
    
    if (!(contrast_name %in% colnames(design)) || !(other_contrast_name %in% colnames(design))) {
      next
    }
    
    contrast <- makeContrasts(contrasts = paste0(contrast_name, " - ", other_contrast_name), 
                              levels = design)
    
    
    fit2 <- eBayes(contrasts.fit(fit, contrast))
    
    sig_genes <- topTable(fit2, number = Inf) %>%
      rownames_to_column("gene") %>%
      filter(adj.P.Val < 0.05 & logFC > log2fc) %>%
      mutate(group1 = celltype,
             group2 = other)
    
    significant_genes <- bind_rows(significant_genes, sig_genes)
  }
  
  return(significant_genes)
}

# Perform comparisons for all cell types
all_comparisons <- data.frame()


for (celltype in celltypes) {
  other_celltypes <- celltypes[celltypes != celltype]
  comparison_results <- compare_celltype(celltype, other_celltypes, design, fit)
  all_comparisons <- bind_rows(all_comparisons, comparison_results)
}

# Display the first few rows of the results
print(head(all_comparisons))
#save
write.csv(all_comparisons, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/all_comparisons.csv")

#I have a table called all_comparisons with columns "gene"      "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"        "group1"    "group2"
#Below I defined different groupings of comparisons. Here I am trying to find the genes that are significant in all comparisons of a grouping, and in the same direction.
#Give me the corresponding gene_lists, and whether they are up or down regulated in each comparison
print(head(coef(fit)))
#DE groupings - MAKE SURE direction of changes is same between members
grouping_1 <- list(list("Stem_cells"), list("Melanotrophs","Corticotrophs",
                                            "Somatotrophs", "Lactotrophs","Thyrotrophs",
                                            "Gonadotrophs"
))

grouping_2 <- list(list("Gonadotrophs"), list("Melanotrophs","Corticotrophs",
                                              "Somatotrophs", "Lactotrophs","Thyrotrophs",
                                              "Stem_cells"
                                              
))

grouping_3 <- list(list("Melanotrophs","Corticotrophs"), list(
  "Somatotrophs", "Lactotrophs","Thyrotrophs",
  "Gonadotrophs", "Stem_cells"
))

grouping_4 <- list(list("Melanotrophs"), list("Corticotrophs"))

grouping_5 <- list(list("Somatotrophs", "Lactotrophs","Thyrotrophs"), list(
  "Gonadotrophs", "Stem_cells", "Melanotrophs",
  "Corticotrophs"))

grouping_6 <- list(list("Lactotrophs"), list("Somatotrophs", "Thyrotrophs"))

grouping_7 <- list(list("Somatotrophs"), list("Lactotrophs", "Thyrotrophs"))

grouping_8 <- list(list("Thyrotrophs"), list("Lactotrophs", "Somatotrophs"))


#############
#Processing groupings
#############
# Function to find significant genes in the same direction across all comparisons in a grouping
find_significant_genes <- function(all_comparisons, groupings1, groupings2, adj_p_threshold = 0.05) {
  # Create all pairwise comparisons
  comps_left <- c()
  comps_right <- c()
  n_comps <- length(groupings1) * length(groupings2)
  
  # Corrected nested loops
  for (i in 1:length(groupings1)) {
    for (j in 1:length(groupings2)) {
      comps_left <- c(comps_left, groupings1[[i]])
      comps_right <- c(comps_right, groupings2[[j]])
    }
  }
  
  # Function to get significant genes for a single comparison
  get_sig_genes <- function(comparison) {
    subset <- all_comparisons[all_comparisons$group1 == comparison$group1 & 
                                all_comparisons$group2 == comparison$group2, ]
    sig_genes <- subset$gene[subset$adj.P.Val < adj_p_threshold]
    
    data.frame(gene = sig_genes, 
               direction = ifelse(subset$logFC[subset$gene %in% sig_genes] > 0, "up", "down"),
               log2fc = subset$logFC[subset$gene %in% sig_genes],
               pvalue = subset$adj.P.Val[subset$gene %in% sig_genes] + 10 ^ -300,
               AveExpr = subset$AveExpr[subset$gene %in% sig_genes],
               stringsAsFactors = FALSE)
  }
  
  # Get significant genes for all comparisons and stitch to a single df
  sig_genes_list <- data.frame(gene = character(), direction = character(), stringsAsFactors = FALSE)
  for (i in 1:length(comps_left)) {
    new_genes <- get_sig_genes(data.frame(group1 = comps_left[i], group2 = comps_right[i]))
    sig_genes_list <- rbind(sig_genes_list, new_genes)
  }
  
  # Only keep genes that occur in all comparisons - meaning n_comps times
  sig_genes_list <- sig_genes_list %>% 
    group_by(gene) %>%
    filter(n() == n_comps) %>%
    ungroup()
  
  # Only keep genes where the direction is the same in all comparisons
  sig_genes_list <- sig_genes_list %>% 
    group_by(gene) %>%
    filter(n_distinct(direction) == 1) %>%
    ungroup()
  
  #add column mean_log2fc
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% mutate(mean_log2fc = mean(log2fc))
  #add geom_mean_adj_pval
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% mutate(geom_mean_adj_pval = exp(mean(log(pvalue))))
  #add mean AvgExpr
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% mutate(mean_AvgExpr = mean(AveExpr))
  #keep only unique rows gene_name and grouping, but keep other rows as well
  sig_genes_list <- sig_genes_list %>% distinct(gene, direction, .keep_all = TRUE)
  
  return(sig_genes_list)
}

# Usage
df <- find_significant_genes(all_comparisons, grouping_5[[1]], grouping_5[[2]])
print(df)


# List of all groupings
all_groupings <- list(grouping_1, grouping_2, grouping_3, grouping_4, 
                      grouping_5, grouping_6, grouping_7, grouping_8)

# Function to process a single grouping
process_grouping <- function(grouping, grouping_name) {
  result <- find_significant_genes(all_comparisons, grouping[[1]], grouping[[2]])
  if (nrow(result) > 0) {
    result$grouping <- grouping_name
    return(result)
  } else {
    return(data.frame(gene = character(0), direction = character(0), grouping = character(0)))
  }
}

# Process all groupings
results_list <- lapply(seq_along(all_groupings), function(i) {
  process_grouping(all_groupings[[i]], paste0("grouping_", i))
})

# Combine all results into a single dataframe
all_results <- do.call(rbind, results_list)

#/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/epitome/data/gene_group_annotation/v_0.01/cpdb.csv
cpdb <- read_csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/gene_group_annotation/v_0.01//cpdb.csv")
#keep where category is TF and keep gene
tf_list <- cpdb %>% filter(category == "TF") %>% select(gene) %>% distinct()
tf_list <- tf_list$gene
options(dplyr.print_max = 50)

for (i in 1:8) {
  grouping_results <- all_results %>% filter(grouping == paste0("grouping_", i))
  cat("\nGrouping", i, "results:\n")
  cat("Number of significant genes:", nrow(grouping_results), "\n")
  if (nrow(grouping_results) > 0) {
    print("Tfs up")
    # Print up to 20 rows by using print() with n parameter
    print(grouping_results %>% 
            filter(gene %in% tf_list & direction == "up") %>% 
            arrange(geom_mean_adj_pval),
          n = 50)
    
    print("Tfs down")
    print(grouping_results %>% 
            filter(gene %in% tf_list & direction == "down") %>% 
            arrange(geom_mean_adj_pval),
          n = 50)
  }
}


#add a column to grouping_results called TF which is 1 if in tf_list
all_results <- all_results %>% mutate(TF = ifelse(gene %in% tf_list, 1, 0))
table(all_results$direction,all_results$grouping)
#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
write_csv(all_results, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/grouping_lineage_markers_0909.csv")
write_csv(all_results, paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/markers/",version,"/grouping_lineage_markers.csv"))
write_csv(all_results, paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/heatmap/",version,"/rna_grouping_lineage_markers.csv"))


#load lineage markers
all_results <- read_csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/grouping_lineage_markers_0909.csv")
table(all_results$grouping,all_results$direction)
        

library(edgeR)
library(Matrix)



# Apply normalization
normalize_expression <- function(expression_table_original, dge_original) {
  norm_factors <- dge_original$samples$norm.factors
  lib_sizes <- dge_original$samples$lib.size
  
  # divide by effective library size
  normalized_expression <- t(t(expression_table_original) / (norm_factors*lib_sizes))
  # Multiply back by 10^6
  normalized_expression <- t(t(normalized_expression) * 10^6)
  
  return(normalized_expression)
}



# Normalize the data
normalized_expression <- normalize_expression(expression_table_original, dge_original)

writeMM(normalized_expression, paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/expression/",version,"/normalized_data.mtx"))
write_csv(meta_data_original, paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/expression/",version,"/meta_data.csv"))
write.table(rownames(normalized_expression), paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/expression/",version,"/genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(normalized_expression), paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/expression/",version,"/datasets.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

rownames(normalized_expression)



##########################################
#make PCA of normalized_expression. Plot PC1 and PC2 twice and color by assignment and exp_type
##########################################

#remove top 20 most highly expressed genes
keep_genes <- rownames(normalized_expression)[order(apply(normalized_expression, 1, sum), decreasing = TRUE)][-(1:20)]
#remove top 10
normalized_expression <- normalized_expression[keep_genes, ]
#keep genes with at least 20000 counts
keep_genes <- rownames(normalized_expression)[apply(normalized_expression, 1, sum) > 20000]
normalized_expression <- normalized_expression[keep_genes, ]
#keep genes with top 2000 highest coefficeints of variation
cv <- apply(normalized_expression, 1, sd) / apply(normalized_expression, 1, mean)
top_genes <- names(sort(cv, decreasing = TRUE)[1:3000])
normalized_expression <- normalized_expression[top_genes, ]

# Load required libraries
library(ggplot2)
library(dplyr)

# Perform PCA
pca_data <- prcomp(t(log1p(normalized_expression)), center = TRUE, scale. = TRUE)

# Calculate variance explained
var_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100

# Create dataframe for plotting
plot_data <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  Assignment = meta_data$assignments,
  ExpType = meta_data$exptype
)

# Create plot by Assignment
p1 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Assignment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA by Cell Assignment"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

# Create plot by ExpType
p2 <- ggplot(plot_data, aes(x = PC1, y = PC2, fill = ExpType)) +
  geom_point(size = 4, alpha = 0.5, shape = 21, stroke = 0.5, color = "black") +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA by Experiment Type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )



# Display plots
print(p1)
print(p2)

#ggsave as png and svg to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/PCA_by_assignment.png", p2, width = 6, height = 6, dpi = 300)
#and svg
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/PCA_by_assignment.svg", p2, width = 6, height = 6)





plot_data <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  Assignment = meta_data$assignments,
  Sex = meta_data$Comp_sex
)

# Create plot by Assignment
p1 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Assignment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA by Sex"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

# Create plot by ExpType
p2 <- ggplot(plot_data, aes(x = PC1, y = PC2, fill = Sex)) +
  geom_point(size = 4, alpha = 0.5, shape = 21, stroke = 0.5, color = "black") +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA by Sex"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )



# Display plots
print(p1)
print(p2)

#ggsave as png and svg to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/PCA_by_sex.png", p2, width = 6, height = 6, dpi = 300)
#and svg
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/PCA_by_sex.svg", p2, width = 6, height = 6)




# Extract first 20 PCs
pcs_for_umap <- pca_data$x[,1:10]

# Run UMAP
#
library(umap)
umap_result <- umap(pcs_for_umap)

# Create plot data
umap_plot_data <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  assignments = meta_data$assignments
)

# Plot
ggplot(umap_plot_data, aes(x = UMAP1, y = UMAP2, color = assignments)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "UMAP Visualization") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )


#repeat pca only on sc subset

normalized_expression_sc <- normalized_expression[, meta_data$exptype == "sc"]
keep_genes_sc <- rowSums(normalized_expression_sc > 1) >= 0.05 * ncol(normalized_expression_sc) & rowSums(normalized_expression_sc) >= 5000
normalized_expression_sc <- normalized_expression_sc[keep_genes_sc, ]


pca_data_sc <- prcomp(t(log1p(normalized_expression_sc)), center = TRUE, scale. = TRUE)
var_explained_sc <- pca_data_sc$sdev^2 / sum(pca_data_sc$sdev^2) * 100
plot_data_sc <- data.frame(
  PC1 = pca_data_sc$x[,1],
  PC2 = pca_data_sc$x[,2],
  Assignment = meta_data$assignments[meta_data$exptype == "sc"]
)

p1 <- ggplot(plot_data_sc, aes(x = PC1, y = PC2, color = Assignment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA by Cell Assignment"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

p1


#install.packages("umap")
library(umap)

# Extract first 20 PCs
pcs_for_umap <- pca_data_sc$x[,1:10]

# Run UMAP
#set 42 seed
set.seed(42)
umap_result <- umap(pcs_for_umap)

# Create plot data
umap_plot_data <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  assignments = meta_data$assignments[meta_data$exptype == "sc"]
)

# Plot
p1 <- ggplot(umap_plot_data, aes(x = UMAP1, y = UMAP2, fill = assignments)) +
  geom_point(size = 4, alpha = 0.5, shape = 21, stroke = 0.5, color = "black") +  # Filled color, thin black outline, transparency
  theme_minimal() +
  labs(title = "UMAP Visualization of pseudobulk samples (sc)") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

p1
#ggsave as png and svg to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/umap_sc.png", p1, width = 6, height = 6, dpi = 300)
#and svg
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/umap_sc.svg", p1, width = 6, height = 6)





normalized_expression_sc <- normalized_expression[, meta_data$exptype == "sn"]
keep_genes_sc <- rowSums(normalized_expression_sc > 1) >= 0.05 * ncol(normalized_expression_sc) & rowSums(normalized_expression_sc) >= 5000
normalized_expression_sc <- normalized_expression_sc[keep_genes_sc, ]


pca_data_sc <- prcomp(t(log1p(normalized_expression_sc)), center = TRUE, scale. = TRUE)
var_explained_sc <- pca_data_sc$sdev^2 / sum(pca_data_sc$sdev^2) * 100
plot_data_sc <- data.frame(
  PC1 = pca_data_sc$x[,1],
  PC2 = pca_data_sc$x[,2],
  Assignment = meta_data$assignments[meta_data$exptype == "sn"]
)

p1 <- ggplot(plot_data_sc, aes(x = PC1, y = PC2, color = Assignment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA by Cell Assignment"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

p1


#install.packages("umap")
library(umap)

# Extract first 20 PCs
pcs_for_umap <- pca_data_sc$x[,1:10]

# Run UMAP
#set 42 seed
set.seed(42)
umap_result <- umap(pcs_for_umap)

# Create plot data
umap_plot_data <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  assignments = meta_data$assignments[meta_data$exptype == "sn"]
)

# Plot
p1 <- ggplot(umap_plot_data, aes(x = UMAP1, y = UMAP2, fill = assignments)) +
  geom_point(size = 4, alpha = 0.5, shape = 21, stroke = 0.5, color = "black") +  # Filled color, thin black outline, transparency
  theme_minimal() +
  labs(title = "UMAP Visualization of pseudobulk samples (sn)") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

p1
#ggsave as png and svg to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/umap_sn.png", p1, width = 6, height = 6, dpi = 300)
#and svg
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/umap_sn.svg", p1, width = 6, height = 6)








normalized_expression_sc <- normalized_expression[, meta_data$exptype == "multi_rna"]
keep_genes_sc <- rowSums(normalized_expression_sc > 1) >= 0.05 * ncol(normalized_expression_sc) & rowSums(normalized_expression_sc) >= 5000
normalized_expression_sc <- normalized_expression_sc[keep_genes_sc, ]


pca_data_sc <- prcomp(t(log1p(normalized_expression_sc)), center = TRUE, scale. = TRUE)
var_explained_sc <- pca_data_sc$sdev^2 / sum(pca_data_sc$sdev^2) * 100
plot_data_sc <- data.frame(
  PC1 = pca_data_sc$x[,1],
  PC2 = pca_data_sc$x[,2],
  Assignment = meta_data$assignments[meta_data$exptype == "multi_rna"]
)

p1 <- ggplot(plot_data_sc, aes(x = PC1, y = PC2, color = Assignment)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA by Cell Assignment"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

p1


#install.packages("umap")
library(umap)

# Extract first 20 PCs
pcs_for_umap <- pca_data_sc$x[,1:10]

# Run UMAP
#set 42 seed
set.seed(42)
umap_result <- umap(pcs_for_umap)

# Create plot data
umap_plot_data <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  assignments = meta_data$assignments[meta_data$exptype == "multi_rna"]
)

# Plot
p1 <- ggplot(umap_plot_data, aes(x = UMAP1, y = UMAP2, fill = assignments)) +
  geom_point(size = 4, alpha = 0.5, shape = 21, stroke = 0.5, color = "black") +  # Filled color, thin black outline, transparency
  theme_minimal() +
  labs(title = "UMAP Visualization of pseudobulk samples (multiome)") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )



p1
#ggsave as png and svg to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/umap_multi.png", p1, width = 6, height = 6, dpi = 300)
#and svg
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/umap_multi.svg", p1, width = 6, height = 6)







##################
##################
##################
#Looking at genes changing with time - age-dependent expression
#####
#####

library(DESeq2)
library(Seurat)
library(tidyverse)
#update R

library(reticulate)
# Specify the path to the .h5ad file
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/pdatas0828.h5ad"

# Read the .h5ad file using reticulate
anndata <- import("anndata")
py_index <- import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)


# Convert the data to Seurat object
# Seurat expects the count matrix to be in column-major order (genes x cells)
counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())
vars

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = counts, meta.data = meta_data)
#print unique Author
rownames(seurat_object) <- vars
seurat_object$assignments <- seurat_object$new_cell_type

#normalize
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1000000)

root = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/"
figs_folder = paste(root,"Figures/",sep="")


meta_data = seurat_object@meta.data
print(length(unique(meta_data$Author)))
print(length(unique(meta_data$SRA_ID)))
#print sum psbulk_n_cells
print(sum(meta_data$psbulk_n_cells))
expression_table = seurat_object@assays$RNA$counts

###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$Modality)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)


#remove where assignment is Blood
expression_table <- expression_table[,meta_data$assignments != "Erythrocytes"]
meta_data <- meta_data[meta_data$assignments != "Erythrocytes",]


meta_data$Age_numeric <- as.numeric(as.character(meta_data$Age_numeric))
meta_data$Comp<- as.numeric(as.character(meta_data$Comp_sex))

#shift to 0 and log10
#remove any with values before 0 
expression_table <- expression_table[,meta_data$Age_numeric >= 0]
meta_data <- meta_data[meta_data$Age_numeric >= 0,]

#keep only where
#meta_data$Age_numeric <- meta_data$Age_numeric - min(meta_data$Age_numeric)
meta_data$Age_numeric <- log10(meta_data$Age_numeric)

#remove where Normal != 1
#expression_table <- expression_table[,meta_data$Normal == 1]
#meta_data <- meta_data[meta_data$Normal == 1,]


# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$Modality)
meta_data$Sex <- make.names(meta_data$Comp)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)
meta_data$Sex <- gsub("[^A-Za-z0-9_]", "", meta_data$Sex)


#add a new column called index
meta_data$index <- 1:nrow(meta_data)
#make index the rownames
rownames(meta_data) <- meta_data$index

#divide values in each column by respective value in meta_data meta_data$psbulk_n_cells
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")

expression_table <- expression_table / meta_data$psbulk_n_cells
#multiply back by 1000
expression_table <- expression_table * median(meta_data$psbulk_n_cells)

#make boxplot of colsums log y
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")


exp_table = expression_table[,meta_data$exptype == "sc"]
meta_data = meta_data[meta_data$exptype == "sc",]

dge = DGEList(counts = exp_table, group = meta_data$assignments)

# Create the design matrix
design <- model.matrix(~0 + assignments +  assignments:Age_numeric, data = meta_data)


#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 3000, min.prop= 0.2)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)

#print genes names with top counts
rownames(dge)[order(rowSums(dge$counts), decreasing = TRUE)[1:15]]
#also print their total count
rowSums(dge$counts)[order(rowSums(dge$counts), decreasing = TRUE)[1:15]]


lib_sizes = colSums(dge$counts[-order(rowSums(dge$counts), decreasing = TRUE)[1:15],])

#add these lib_size to dge
dge$lib.size = lib_sizes

#norm
dge <- calcNormFactors(dge, method = "TMM")

colnames(design)<-make.names(colnames(design))

fit <- voom(dge, design, plot=TRUE)
fit <- lmFit(fit, design)

#####
#Age effect
#####

contrast <- makeContrasts(
  #assignmentsStem_cells - assignmentsLactotrophs, 
  assignmentsStem_cells.Age_numeric,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


cpdb <- read_csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/gene_group_annotation/v_0.01//cpdb.csv")
#keep where category is TF and keep gene
tf_list <- cpdb %>% filter(category == "TF") %>% select(gene) %>% distinct()
tf_list <- tf_list$gene


print("Up TFs with time")
#rank it by adj.P.val
print(top_genes %>% filter(genes %in% tf_list & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))
print("Down TFs with time")
print(top_genes %>% filter(genes %in% tf_list & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))


print("Up with time")
#rank it by adj.P.val
print(top_genes %>% filter( logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))
print("Down  with time")
print(top_genes %>% filter( logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))



#save these tfs up and down
tfs_up = top_genes %>% filter(genes %in% tf_list & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
tfs_down = top_genes %>% filter(genes %in% tf_list & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)

tfs_up

tfs_down

#save to "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
write.csv(tfs_up, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_tfs_up.csv")
write.csv(tfs_down, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_tfs_down.csv")

cpdb <- read_csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/gene_group_annotation/v_0.01//cpdb.csv")
#keep where category is TF and keep gene
ligands  <- cpdb %>% filter(category == "ligand") %>% select(gene) %>% distinct()
ligands <- ligands$gene
receptors  <- cpdb %>% filter(category == "receptor") %>% select(gene) %>% distinct()
receptors <- receptors$gene

#ligands /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/CellChatDBmouse.csv
ligands_up = top_genes %>% filter(genes %in% ligands & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
ligands_down = top_genes %>% filter(genes %in% ligands & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
ligands_up
ligands_down


#receptors
receptors_up = top_genes %>% filter(genes %in% receptors & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
receptors_down = top_genes %>% filter(genes %in% receptors & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
receptors_up
receptors_down

write.csv(ligands_up, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_ligands_up.csv")
write.csv(ligands_down, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_ligands_down.csv")
write.csv(receptors_up, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_receptors_up.csv")
write.csv(receptors_down, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_receptors_down.csv")


library("clipr")
up_genes <- top_genes[top_genes$logFC > 0 & top_genes$adj.P.Val < 0.05,]
#print top 10 highest pval


print(head(up_genes,10))
write_clip(paste(up_genes$gene, collapse = "\n"))
down_genes <- top_genes[top_genes$logFC < 0 & top_genes$adj.P.Val < 0.05,]
#print top 10 highest pval
print(head(down_genes,10))
write_clip(paste(down_genes$gene, collapse = "\n"))




#Extracting for all cell types



coefs = c("assignmentsCorticotrophs.Age_numeric",
          "assignmentsEndothelial_cells.Age_numeric", "assignmentsGonadotrophs.Age_numeric",     
          "assignmentsImmune_cells.Age_numeric",      "assignmentsLactotrophs.Age_numeric",      
          "assignmentsMelanotrophs.Age_numeric",      "assignmentsMesenchymal_cells.Age_numeric",
          "assignmentsPituicytes.Age_numeric",        "assignmentsSomatotrophs.Age_numeric",     
          "assignmentsStem_cells.Age_numeric",        "assignmentsThyrotrophs.Age_numeric" )


# Create empty list to store results
all_results <- list()

# Extract cell types from coefficient names
get_cell_type <- function(coef_name) {
  # Remove "assignments" prefix and ".Age_numeric" suffix
  gsub("assignments(.*)\\.Age_numeric", "\\1", coef_name)
}

# Iterate through each coefficient
for(coef in coefs) {
  # Create contrast matrix
  contrast <- makeContrasts(
    contrasts = coef,
    levels = design
  )
  
  # Compute differential expression
  diff_expression <- contrasts.fit(fit, contrast)
  diff_expression <- eBayes(diff_expression, robust=TRUE)
  
  # Get top genes
  top_genes <- topTable(diff_expression, number = Inf)
  top_genes$genes <- rownames(top_genes)
  
  # Add cell type column
  top_genes$cell_type <- get_cell_type(coef)
  
  # Store results
  all_results[[coef]] <- top_genes
}

# Combine all results into single dataframe
final_results <- do.call(rbind, all_results)

# Reset row names
rownames(final_results) <- NULL



#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/aging/v_0.01
version = "v_0.01"
write.csv(final_results, paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/aging/v_0.01/aging_genes.csv"))


final_results <- read.csv(paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/aging/v_0.01/aging_genes.csv"), header = TRUE, stringsAsFactors = FALSE)

#keep only padj < 0.05
final_results_pan_pit <- final_results[final_results$adj.P.Val < 0.05,]
#remove where cell_type "Mesenchymal_cells" "Pituicytes"   "Immune_cells" "Endothelial_cells"
final_results_pan_pit <- final_results_pan_pit[!final_results_pan_pit$cell_type %in% c("Mesenchymal_cells", "Pituicytes", "Immune_cells", "Endothelial_cells"),]
#calc avg logFC and geom avg adj pval
final_results_pan_pit$avg_logFC <- ave(final_results_pan_pit$logFC, final_results_pan_pit$genes, FUN = mean)
#geom avg adj pval
final_results_pan_pit$avg_adj_pval <- ave(final_results_pan_pit$adj.P.Val, final_results_pan_pit$genes, FUN = function(x) prod(x)^(1/length(x)))
#add another column saying number of times a given gene occurs
final_results_pan_pit$gene_count <- ave(final_results_pan_pit$genes, final_results_pan_pit$genes, FUN = length)

final_results_pan_pit <- final_results_pan_pit[!duplicated(final_results_pan_pit$genes),]
#at least 3 gene count
final_results_pan_pit <- final_results_pan_pit[final_results_pan_pit$gene_count >= 3,]
final_results_pan_pit

#order according to gene_count
final_results_pan_pit <- final_results_pan_pit[order(final_results_pan_pit$gene_count, decreasing = TRUE),]


#I have this df called final_results with columns logFC, genes, adj.P.Val and cell_type. Make a heatmap, where tiles are outline wiht red if significant and colored with light to dark blue based on logFC. Make columns cell_types and rows genes. Specifically look at Ifit3, Mki67 and Lef1


# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)

# Highlight specific genes of interest
genes_of_interest <- c("Wnt10a","Fzd2","Fzd9","Wif1","Lef1","Rspo4","Sfrp1","Sfrp2","Sfrp5","Apcdd1","Tnfrsf19","Sostdc1",
                       "Cdk1","H19","E2f8","Mki67","Top2a","Ube2c","Aurkb",
                       "Lcn2","Cxcl13","Ccl5","H2-Aa","H2-Ab1","Ifit3","Ifit3b","C1qc","Il1b","Il6")
                       
                       
# Create the heatmap
plot <- final_results %>%
  # Filter to only include genes of interest
  filter(genes %in% genes_of_interest) %>%
  # Create a significance flag and preserve input order
  mutate(is_significant = adj.P.Val < 0.05,
         # Add line thickness variable
         line_size = ifelse(genes %in% genes_of_interest, 0.6, 0.4),
         # Convert genes to factor with levels in the exact order specified
         genes = factor(genes, levels = genes_of_interest)) %>%
  # Create the plot
  ggplot(aes(x = cell_type, y = genes, fill = logFC)) +
  geom_tile(aes(color = is_significant, size = line_size), 
            width = 0.8, height = 0.8) +
  # Set color scales
  scale_fill_gradient2(low = "#ffa500", mid = "white", 
                       high = "#63b3ed", midpoint = 0) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                     name = "Significant",
                     labels = c("FALSE" = "Not Significant", "TRUE" = "Significant"),
                     guide = guide_legend(override.aes = list(
                       fill = "white",
                       size = 1,
                       linetype = 1,
                       linewidth = 3
                     ))) +
  scale_size_identity() +
  # Customize theme
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  #y axis italics
  theme(axis.text.y = element_text(size = 10, face = "italic")) +
  labs(fill = "logFC")

plot


#save plot to/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/DE
ggsave(plot, filename = "/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/DE/heatmap_aging_genes.png", width = 5.1, height = 5.0, dpi = 300)
#also svg
ggsave(plot, filename = "/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/DE/heatmap_aging_genes.svg", width = 5.1, height = 5.0, dpi = 300)



# Process the final_results data for plotting
library(ggplot2)
library(dplyr)

# Create a dataframe to store the counts of age-related genes per cell type
age_genes <- list()
cell_types <- unique(final_results$cell_type)

for (cell_type in cell_types) {
  # Filter results for this cell type
  cell_results <- final_results[final_results$cell_type == cell_type, ]
  
  # Count significant genes (adj.P.Val < 0.05)
  # Positive logFC = increased with age, negative = decreased with age
  increased_genes <- nrow(cell_results[cell_results$adj.P.Val < 0.05 & cell_results$logFC > 0, ])
  decreased_genes <- nrow(cell_results[cell_results$adj.P.Val < 0.05 & cell_results$logFC < 0, ])
  
  age_genes[[cell_type]] <- list(
    increased = increased_genes,
    decreased = decreased_genes
  )
  
  cat(sprintf("%s: %d increased, %d decreased genes with aging\n", 
              cell_type, increased_genes, decreased_genes))
}

# Create dataframe for plotting
plot_data <- data.frame(
  cell_type = names(age_genes),
  increased_genes = sapply(age_genes, function(x) x$increased),
  decreased_genes = sapply(age_genes, function(x) x$decreased)
)
#remove Immune_cells
plot_data <- plot_data[!grepl("Immune_cells", plot_data$cell_type), ]
# Reshape the data for plotting
plot_data_long <- rbind(
  data.frame(
    cell_type = plot_data$cell_type,
    count = -plot_data$decreased_genes,  # Negative for decreased to show on left
    change = "Decreased",
    label = plot_data$decreased_genes    # Positive label
  ),
  data.frame(
    cell_type = plot_data$cell_type,
    count = plot_data$increased_genes,
    change = "Increased",
    label = plot_data$increased_genes
  )
)

# Format labels
plot_data_long$label[plot_data_long$label == 0] <- ""
plot_data_long$label <- ifelse(plot_data_long$change == "Increased", 
                               paste0("   ", plot_data_long$label),
                               plot_data_long$label)

# Ensure count is numeric
plot_data_long$count <- as.numeric(plot_data_long$count)

# Calculate total differentially expressed genes per cell type and reorder
plot_data_long <- plot_data_long %>%
  group_by(cell_type) %>%
  mutate(total_diff_exp = sum(abs(count))) %>% # Total absolute differentially expressed genes
  ungroup() %>%
  arrange(desc(total_diff_exp)) %>%            # Order by total genes
  mutate(cell_type = factor(cell_type, levels = unique(cell_type))) # Maintain order

# Add a column for fixed label positions
plot_data_long <- plot_data_long %>%
  mutate(label_position = ifelse(count > 0, 400, -400)) # Labels at fixed positions

# Plot
plt = ggplot(plot_data_long, aes(x = cell_type, y = count, fill = change)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = label, y = label_position), # Use fixed label positions
            color = "black",
            size = 5) +
  scale_fill_manual(values = c("Decreased" = "#FFA500", "Increased" = "#63B3ED")) +
  scale_y_continuous(
    labels = abs,  # Show absolute values on axis
    limits = c(-max(abs(plot_data_long$count)) * 1.2, max(abs(plot_data_long$count)) * 1.2) # Symmetric axis
  ) +
  coord_flip() +
  labs(
    title = "Age-related Genes Across Cell Types",
    x = "Cell Type",
    y = "Number of Age-related Genes",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16)
  )

print(plt)




# Save the plot as PNG and SVG
output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE"

# PNG with high resolution
ggsave(
  filename = file.path(output_dir, "age_related_genes.png"),
  plot = plt,
  width = 6,
  height = 5,
  dpi = 300
)

# SVG format
ggsave(
  filename = file.path(output_dir, "age_related_genes.svg"),
  plot = plt,
  width = 6,
  height = 5
)

print(plt)






#############
############
#############
############
#############

##################
##################
##################
#Looking at genes changing with sex

library(DESeq2)
library(Seurat)
library(tidyverse)
#update R
library(reticulate)
# Specify the path to the .h5ad file
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/pdatas0828.h5ad"

# Read the .h5ad file using reticulate
anndata <- reticulate::import("anndata")
py_index <- reticulate::import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)


# Convert the data to Seurat object
# Seurat expects the count matrix to be in column-major order (genes x cells)
counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())
vars

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = counts, meta.data = meta_data)
#print unique Author
rownames(seurat_object) <- vars
seurat_object$assignments <- seurat_object$new_cell_type

#normalize
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1000000)


root = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/"
figs_folder = paste(root,"Figures/",sep="")


meta_data = seurat_object@meta.data
print(length(unique(meta_data$Author)))
print(length(unique(meta_data$SRA_ID)))
#print sum psbulk_n_cells
print(sum(meta_data$psbulk_n_cells))
expression_table = seurat_object@assays$RNA$counts

###################################################
### Now with Limma-voom below
library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$Modality)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)


#remove where assignment is Blood
expression_table <- expression_table[,meta_data$assignments != "Erythrocytes"]
meta_data <- meta_data[meta_data$assignments != "Erythrocytes",]


meta_data$Age_numeric <- as.numeric(as.character(meta_data$Age_numeric))
meta_data$Comp<- as.numeric(as.character(meta_data$Comp_sex))
#shift to 0 and log10
#remove any with values before 0 
expression_table <- expression_table[,meta_data$Age_numeric >= 10]
meta_data <- meta_data[meta_data$Age_numeric >= 10,]

#keep only samples between 20 and 150
expression_table <- expression_table[,meta_data$Age_numeric <= 150]
meta_data <- meta_data[meta_data$Age_numeric <= 150,]
expression_table <- expression_table[,meta_data$Age_numeric > 10]
meta_data <- meta_data[meta_data$Age_numeric > 10,]

#keep only Normal == 1
#expression_table <- expression_table[,meta_data$Normal == 1]
#meta_data <- meta_data[meta_data$Normal == 1,]

meta_data$exptype[meta_data$`10X version` == "PARSE_WT_MegaV2"] <- "parse"
meta_data$exptype[meta_data$`10X version` == "Parse_WT"] <- "parse"

#remove those that are parse
expression_table <- expression_table[,meta_data$exptype != "parse"]
meta_data <- meta_data[meta_data$exptype != "parse",]


# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$exptype)
meta_data$Sex <- make.names(meta_data$Comp)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)
meta_data$Sex <- gsub("[^A-Za-z0-9_]", "", meta_data$Sex)


#order meta_data based on exptype
#reorder expression_table based on meta_data index
expression_table <- expression_table[,order(meta_data$exptype)]
meta_data <- meta_data[order(meta_data$exptype),]
#print all values
print(meta_data$exptype)
#add a new column called index
meta_data$index <- 1:nrow(meta_data)
#make index the rownames
rownames(meta_data) <- meta_data$index


#divide values in each column by respective value in meta_data meta_data$psbulk_n_cells
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")

expression_table <- expression_table / meta_data$psbulk_n_cells
#multiply back by 1000
expression_table <- expression_table * median(meta_data$psbulk_n_cells)

#make boxplot of colsums log y
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")


exp_table_sc = expression_table[,meta_data$exptype == "sc"]
meta_data_sc = meta_data[meta_data$exptype == "sc",]

exp_table_sn = expression_table[,meta_data$exptype == "sn"]
meta_data_sn = meta_data[meta_data$exptype == "sn",]

exp_table_multi = expression_table[,meta_data$exptype == "multi_rna"]
meta_data_multi = meta_data[meta_data$exptype == "multi_rna",]

#exp_table_parse = expression_table[,meta_data$exptype == "parse"]
#meta_data_parse = meta_data[meta_data$exptype == "parse",]


dge = DGEList(counts = expression_table, group = meta_data$assignments)
# Prepare the DGEList object for limma voom
dge_sc <- DGEList(counts = exp_table_sc, group = meta_data_sc$assignments)
dge_sn <- DGEList(counts = exp_table_sn, group = meta_data_sn$assignments)
dge_multi <- DGEList(counts = exp_table_multi, group = meta_data_multi$assignments)
#dge_parse <- DGEList(counts = exp_table_parse, group = meta_data_parse$assignments)


library(limma)
library(edgeR)
design <- model.matrix(~0 + assignments +  assignments:Sex, data = meta_data)
#make ecdf of total counts of genes
plot(ecdf(rowSums(dge$counts)), main = "ecdf of total counts of genes", xlab = "Total counts", ylab = "Proportion of genes")
#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 3000, min.prop= 0.2)

dge <- dge[keep_genes, ]
dim(dge)
gene_names <- rownames(dge)
dge_sc <- dge_sc[gene_names,]
dge_sn <- dge_sn[gene_names,]
dge_multi <- dge_multi[gene_names,]
#dge_parse <- dge_parse[gene_names,]

lib_sizes_sc = colSums(dge_sc$counts[-order(rowSums(dge_sc$counts), decreasing = TRUE)[1:15],])
lib_sizes_sn = colSums(dge_sn$counts[-order(rowSums(dge_sn$counts), decreasing = TRUE)[1:15],])
lib_sizes_multi = colSums(dge_multi$counts[-order(rowSums(dge_multi$counts), decreasing = TRUE)[1:15],])
#lib_sizes_parse = colSums(dge_parse$counts[-order(rowSums(dge_parse$counts), decreasing = TRUE)[1:15],])

#add these lib_size to dge
dge_sc$samples$lib.size= lib_sizes_sc
dge_sn$samples$lib.size = lib_sizes_sn
dge_multi$samples$lib.size = lib_sizes_multi
#dge_parse$samples$lib.size = lib_sizes_parse

#norm
dge_sc <- calcNormFactors(dge_sc, method = "TMM")
dge_sn <- calcNormFactors(dge_sn, method = "TMM")
dge_multi <- calcNormFactors(dge_multi, method = "TMM")
#dge_parse <- calcNormFactors(dge_parse, method = "TMM")


design_sc <- model.matrix(~ 0 + assignments + assignments:Sex, data = meta_data_sc)
design_sn <- model.matrix(~ 0 + assignments + assignments:Sex, data = meta_data_sn)
design_multi <- model.matrix(~ 0 + assignments + assignments:Sex, data = meta_data_multi)
#design_parse <- model.matrix(~ 0 + assignments + assignments:Sex, data = meta_data_parse)
  

dge$samples = rbind(dge_multi$samples,#dge_parse$samples,
                    dge_sc$samples,dge_sn$samples)


meta_data = rbind(meta_data_multi, #meta_data_parse,
                  meta_data_sc,meta_data_sn)

# Create the fixed effects design matrix
meta_data$assignments_exptype <- paste(meta_data$assignments, meta_data$exptype, sep = "_")

# First convert assignments to factor if it isn't already
meta_data$assignments <- as.factor(meta_data$assignments)

# Create the design matrix
design <- model.matrix(~0 + assignments +  assignments:Sex, data = meta_data)
colnames(design)
#  + assignments:Sex, data = meta_data)
colnames(design)<-make.names(colnames(design))


fit_sc <- voom(dge_sc, design_sc, plot=TRUE)
fit_sn <- voom(dge_sn, design_sn, plot=TRUE)
fit_multi <- voom(dge_multi, design_multi, plot=TRUE)
#fit_parse <- voom(dge_parse, design_parse, plot=TRUE)
fit2=cbind(fit_multi,
           #fit_parse,
           fit_sc,fit_sn)

#double running this based on
#https://support.bioconductor.org/p/114663/


fit <- voom(dge, design, plot=TRUE)
fit$weights <- fit2$weights
corfit <- duplicateCorrelation(fit, design, block=meta_data$exptype)
corfit$consensus.correlation
fit <- voom(dge, design, plot=TRUE, block=meta_data$exptype, correlation=corfit$consensus.correlation)
corfit <- duplicateCorrelation(fit, design, block=meta_data$exptype)
corfit$consensus.correlation

#fit$weights <- matrix(1, nrow=nrow(fit$weights), ncol=ncol(fit$weights))
fit <- lmFit(fit, design, block=meta_data$exptype, correlation=corfit$consensus.correlation)

print(head(coef(fit)))
###
# Make column names valid
colnames(design) <- make.names(colnames(design), unique = TRUE)

print(head(coef(fit)))
#turn coef(fit) into a dataframe
coef_df = as.data.frame(coef(fit))

coef_df$gene = rownames(coef_df)

#####
#Sex effect
#####
contrast <- makeContrasts(
  assignmentsStem_cells.SexX1,
  levels = design
)


# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))

cpdb <- read_csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/gene_group_annotation/v_0.01//cpdb.csv")
#keep where category is TF and keep gene
tf_list <- cpdb %>% filter(category == "TF") %>% select(gene) %>% distinct()
tf_list <- tf_list$gene
tf_list

print("Up TFs with male")
#rank it by adj.P.val
print(top_genes %>% filter(genes %in% tf_list & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))
print("Down TFs with female")
print(top_genes %>% filter(genes %in% tf_list & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))

library(clipr)
# Write the genes to clipboard, one per line
up_genes <- top_genes %>% filter(logFC > 0 & adj.P.Val < 0.05)
clipr::write_clip(paste(up_genes$gene, collapse = "\n"))
down_genes <- top_genes %>% filter(logFC < 0 & adj.P.Val < 0.05)
clipr::write_clip(paste(down_genes$gene, collapse = "\n"))
background_genes <- rownames(dge)
clipr::write_clip(paste(background_genes, collapse = "\n"))


#save these to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_up_genes.csv
write.csv(up_genes, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_sexually_dimorphic_up_genes.csv")
write.csv(down_genes, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_sexually_dimorphic_down_genes.csv")



# Initialize a list to store results
sex_specific_genes <- list()
cell_types <- unique(meta_data$assignments)

for (cell_type in cell_types) {
  # Create contrast for this cell type's sex effect
  contrast <- makeContrasts(
    sex_effect = paste0("assignments", cell_type, ".SexX1"),
    levels = design
  )
  
  # Compute differential expression
  diff_expr <- contrasts.fit(fit, contrast)
  diff_expr <- eBayes(diff_expr, robust=TRUE)
  
  # Get results
  results <- topTable(diff_expr, number=Inf)
  
  #filter out this where Avg_exp is 0
  results <- results[results$AveExpr > 0,]
  
  # Store significant genes (positive logFC = male, negative = female)
  male_genes <- rownames(results[results$adj.P.Val < 0.05 & results$logFC > 0,])
  female_genes <- rownames(results[results$adj.P.Val < 0.05 & results$logFC < 0,])
  
  sex_specific_genes[[cell_type]] <- list(
    male = male_genes,
    female = female_genes
  )
  
  cat(sprintf("%s: %d male-specific, %d female-specific genes\n", 
              cell_type, length(male_genes), length(female_genes)))
  
}

# Create dataframe for plotting
plot_data <- data.frame(
  cell_type = names(sex_specific_genes),
  male_genes = sapply(sex_specific_genes, function(x) length(x$male)),
  female_genes = sapply(sex_specific_genes, function(x) length(x$female))
)



# Reshape the data for plotting
plot_data_long <- rbind(
  data.frame(
    cell_type = plot_data$cell_type,
    count = -plot_data$female_genes,  # Negative for female to show on left
    sex = "Female",
    label = plot_data$female_genes    # Positive label
  ),
  data.frame(
    cell_type = plot_data$cell_type,
    count = plot_data$male_genes,
    sex = "Male",
    label = plot_data$male_genes
  )
)

# Create the plot
#where the fill is =, make it ""
plot_data_long$label[plot_data_long$label == 0] <- ""
plot_data_long$label <- ifelse(plot_data_long$sex == "Male", 
                               paste0("   ", plot_data_long$label),
                               plot_data_long$label)

#remove Posterior_pit
plot_data_long <- plot_data_long[plot_data_long$cell_type != "Posterior_pit",]
plot_data_long$count <- as.numeric(plot_data_long$count)
# Ensure count is numeric
plot_data_long$count <- as.numeric(plot_data_long$count)

# Precompute vjust for each row
plot_data_long$vjust <- ifelse(plot_data_long$count > 0, 0.6, -0.6)
library(dplyr)

# Ensure count is numeric
plot_data_long$count <- as.numeric(plot_data_long$count)
library(dplyr)

# Calculate total differentially expressed genes per cell type and reorder
plot_data_long <- plot_data_long %>%
  group_by(cell_type) %>%
  mutate(total_diff_exp = sum(abs(count))) %>% # Total absolute differentially expressed genes
  ungroup() %>%
  arrange(desc(total_diff_exp)) %>%            # Order by total genes
  mutate(cell_type = factor(cell_type, levels = unique(cell_type))) # Maintain order

# Add a column for fixed label positions
plot_data_long <- plot_data_long %>%
  mutate(label_position = ifelse(count > 0, 400, -400)) # Labels at fixed positions

# Plot
plt = ggplot(plot_data_long, aes(x = cell_type, y = count, fill = sex)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = label, y = label_position), # Use fixed label positions
            color = "black",
            size = 5) +
  scale_fill_manual(values = c("Female" = "#FFA500", "Male" = "#63B3ED")) +
  scale_y_continuous(
    labels = abs,  # Show absolute values on axis
    limits = c(-max(abs(plot_data_long$count)) * 1.2, max(abs(plot_data_long$count)) * 1.2) # Symmetric axis
  ) +
  coord_flip() +
  labs(
    title = "Sex-biased Genes Across Cell Types",
    x = "Cell Type",
    y = "Number of Sex-biased Genes",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16),
    #y tick smaller
    axis.text.y = element_text(size = 14)
  )

plt

stem_sex_spec = sex_specific_genes$Stem_cells$female


#save plt as png and svg to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1
fold = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/"
ggsave(plt, filename = paste0(fold,"sex_biased_genes.png"), width = 6, height = 5.5, dpi = 300)
ggsave(plt, filename = paste0(fold,"sex_biased_genes.svg"), width = 6, height = 5.5, dpi = 300)

plt
#save all sex biased gene hits
print(colnames(coef(fit)))
coefs = c("assignmentsCorticotrophs.SexX1", "assignmentsEndothelial_cells.SexX1","assignmentsGonadotrophs.SexX1",
"assignmentsImmune_cells.SexX1", "assignmentsLactotrophs.SexX1", "assignmentsMelanotrophs.SexX1",
"assignmentsMesenchymal_cells.SexX1", "assignmentsPituicytes.SexX1", "assignmentsSomatotrophs.SexX1",
"assignmentsStem_cells.SexX1", "assignmentsThyrotrophs.SexX1")

df <- data.frame()
df_raw <- data.frame()
cell_types <- unique(meta_data$assignments)

for (cell_type in cell_types) {
  # Create contrast for this cell type's sex effect
  contrast <- makeContrasts(
    sex_effect = paste0("assignments", cell_type, ".SexX1"),
    levels = design
  )
  
  # Compute differential expression
  diff_expr <- contrasts.fit(fit, contrast)
  diff_expr <- eBayes(diff_expr, robust=TRUE)
  
  # Get results
  results <- topTable(diff_expr, number=Inf)
  top_genes <- results
  top_genes$cell_type <- gsub("assignments", "", gsub(".SexX1", "", cell_type))
  top_genes$gene <- rownames(top_genes)
  
  #for those where log2fc <0 say female
  top_genes$sex <- ifelse(top_genes$logFC > 0, "male","female")
  top_genes_filtered <- top_genes[top_genes$adj.P.Val < 0.05,]
  
  df = rbind(df, top_genes_filtered)
  df_raw = rbind(df_raw, top_genes)
  
}
  
#add to df a column called "occurs" which is the number of cell types the gene is differentially expressed in
df <- df %>%
  group_by(gene) %>%
  summarize(occurs = n_distinct(cell_type)) %>%
  left_join(df, by = "gene")

df

#write
write.csv(df,"/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes.csv")
write.csv(df,"/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/sex_dimorphism/v_0.01/sexually_dimorphic_genes.csv")

write.csv(df_raw,"/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes_raw.csv")



df_common_female <- df %>%
  filter(sex == "female") %>%
  group_by(gene) %>%
  summarize(n_cell_types = n_distinct(cell_type)) %>%
  filter(n_cell_types >= 3)

# Same for males
df_common_male <- df %>%
  filter(sex == "male") %>%
  group_by(gene) %>%
  summarize(n_cell_types = n_distinct(cell_type)) %>%
  filter(n_cell_types >= 3)


#####
##### Making heatmap of shared sex genes
#####


df_sex_genes_filtered <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes.csv")
#keep only where occurs >=7
df_sex_genes_filtered <- df_sex_genes_filtered %>% filter(occurs > 5)
df_sex_genes <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes_raw.csv")
#keep genes in df_sex_filtered
df_sex_genes <- df_sex_genes %>% filter(gene %in% df_sex_genes_filtered$gene)

df_sex_genes <- df_sex_genes %>%
  group_by(gene) %>%
  mutate(avg_log2fc = mean(logFC, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(avg_log2fc))

#first 20 gene names
genes_of_interest <- unique(df_sex_genes$gene)

final_results <- df_sex_genes

# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)

#-1 the logFC
final_results$logFC <- -final_results$logFC
# Create the heatmap
plot <- final_results %>%
  # Filter to only include genes of interest
  dplyr::filter(gene %in% genes_of_interest) %>%
  # Create a significance flag and preserve input order
  mutate(is_significant = adj.P.Val < 0.05,
         # Add line thickness variable
         line_size = ifelse(gene %in% genes_of_interest, 0.6, 0.4),
         # Convert genes to factor with levels in the exact order specified
         gene = factor(gene, levels = genes_of_interest)) %>%
  # Create the plot
  ggplot(aes(x = cell_type, y = gene, fill = logFC)) +
  geom_tile(aes(color = is_significant, size = line_size), 
            width = 0.8, height = 0.8) +
  
  # Set color scales
  scale_fill_gradient2(low = "#63b3ed", mid = "white", 
                       high = "#ffa500", midpoint = 0) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                     name = "Significant",
                     labels = c("FALSE" = "Not Significant", "TRUE" = "Significant"),
                     guide = guide_legend(override.aes = list(
                       fill = "white",
                       size = 1,
                       linetype = 1,
                       linewidth = 3
                     ))) +
  scale_size_identity() +
  # Customize theme
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1,
                                   vjust = 1, size = 12),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  #y axis italics
  theme(axis.text.y = element_text(size = 10, face = "italic")) +
  labs(fill = "logFC")

plot

#ggsave to /Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes_heatmap.png", 
       plot = plot, 
       width = 6.5, height = 5.5, dpi = 300)

#svg
ggsave("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes_heatmap.svg", 
       plot = plot, 
       width = 6.5, height = 5.5, dpi = 300)




##################
##################
##################
# Cell typing marker
##################
##################
##################


# Main loop
all_markers <- list()
#initialise empty df 
all_markers_df <- data.frame()
celltypes <- unique(meta_data$assignments)
for (celltype in celltypes) {
  other_celltypes <- setdiff(celltypes, celltype)
  significant_genes <- compare_celltype(celltype, other_celltypes, design, fit, 0)
  print(paste("Cell type:", celltype, "Significant genes found:", nrow(significant_genes)))
  # Find genes significant in all but at most one comparison
  min_required <- length(other_celltypes)
  markers <- significant_genes %>% group_by(gene) %>% dplyr::filter(n() >= min_required) 
  print(paste("Markers after filtering for consistency in", celltype, ":", nrow(markers)))
  #only keep if log2fc is ALWAYS > 2
  markers <- markers %>% group_by(gene) %>% dplyr::filter(all(logFC > 2)) %>% dplyr::summarise(log2fc = mean(logFC), pval = exp(mean(log(adj.P.Val))), avg_expr = mean(AveExpr))
  #order by pval
  print(paste("Markers after filtering for log2fc > 2 in", celltype, ":", nrow(markers)))
  markers <- markers %>% dplyr::arrange(pval)
  #add col about cell types
  markers$celltype <- celltype
  
  
  #add markers to all_markers_df
  all_markers_df <- rbind(all_markers_df, markers)
  
  #remove dups
  markers <- unique(markers%>% pull(gene))
  all_markers[[celltype]] <- markers
}


all_markers_df <- all_markers_df %>% dplyr::arrange(pval)


#save it to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/cell_typing_markers.csv
write.csv(all_markers_df, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/cell_typing_markers.csv")

write.csv(all_markers_df, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/markers/v_0.01/cell_typing_markers.csv")

all_markers_df <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/cell_typing_markers.csv")

#print top 5 marker for each cell type 
all_markers_df <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/cell_typing_markers.csv")
coefs_table <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/coef.csv")
rownames(coefs_table) <- coefs_table$X
#remove first col
coefs_table <- coefs_table[, -1]
#for each col name remove assignments in coefs_table
colnames(coefs_table) <- gsub("assignments", "", colnames(coefs_table))
#add a new col to all_markers_df called avg_exp, get this value by finding the cell_type as a column in coefs_table and then accessing the gene (row_index) and get value
all_markers_df$avg_exp <- NA
for (i in 1:nrow(all_markers_df)) {
  gene <- all_markers_df$gene[i]
  celltype <- all_markers_df$celltype[i]
  all_markers_df$avg_exp[i] <- coefs_table[gene, celltype]
}

#keep those with at least 5 > avg_exp
#all_markers_df <- all_markers_df %>% filter(avg_exp > 5)



celltypes <- unique(all_markers_df$celltype)

for (ct in celltypes) {
  print(all_markers_df %>% dplyr::filter(celltype == ct) %>% head(5))
}


# Create final markers dataframe, handling empty marker lists
markers_df <- all_markers_df %>% dplyr::select(gene, celltype, avg_exp, pval) %>% arrange(celltype, pval)

# Remove rows with NA genes
markers_df <- markers_df %>% filter(!is.na(gene))

# Print summary
print(table(markers_df$celltype))


# Print the first few rows of the final markers dataframe
print("First few rows of the markers dataframe:")
print(head(markers_df))

#print top
# Print structure of meta_data and design matrix for debugging
print("Structure of meta_data:")
print(str(meta_data))

print("First few rows of the design matrix:")
print(head(design))


markers<- markers_df %>%
  filter(celltype == "Stem_cells") %>%
  select(gene)
#convert to ordinary list
markers <- markers$gene
markers <- as.character(markers)

#print top5 for each cell type
for (ct in celltypes) {
  print(markers_df %>% filter(celltype == ct) %>% head(5))
}




#############
###Novel markers other than those used in CellAssign
##########

#load /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/spatial panel construction/cellassign_markers_updated.xlsx
library(readxl)
cellassign_markers = read_excel("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/spatial panel construction/cellassign_markers_updated.xlsx")
cellassign_markers = cellassign_markers$Gene
print(dim(markers_df))
new_markers_not_in_cellassign <- markers_df[!markers_df$gene %in% cellassign_markers,]
print(dim(new_markers_not_in_cellassign))

new_markers_not_in_cellassign
#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/novel_markers.csv
file_to_save = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/novel_markers.csv"
write.csv(new_markers_not_in_cellassign, file_to_save)

new_markers_not_in_cellassign <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/novel_markers.csv")
#also remove those with Rik in the name
#new_markers_not_in_cellassign <- new_markers_not_in_cellassign %>% filter(!grepl("Rik", gene))
#remove those with gm or ENSMUSG in the name
#new_markers_not_in_cellassign <- new_markers_not_in_cellassign %>% filter(!grepl("Gm", gene))
#all_markers_df <- all_markers_df %>% filter(!grepl("ENSMUSG", gene))

#avg_exp above 5
new_markers_not_in_cellassign <- new_markers_not_in_cellassign %>% filter(avg_exp > 3)

for (ct in celltypes) {
  print(new_markers_not_in_cellassign %>% filter(celltype == ct) %>% head(6))
}



#keep only where gene is in new_markers_not_in_cellassign
plot_markers_df <- all_markers_df %>% filter(gene %in% new_markers_not_in_cellassign$gene)
#order by lowest pval first
plot_markers_df <- plot_markers_df %>% arrange(pval)




#remove Nr5a1os
plot_markers_df <- plot_markers_df %>% filter(gene != "Nr5a1os")
for (ct in celltypes) {
  print(plot_markers_df  %>% filter(celltype == ct) %>% select(gene,celltype) %>% head(5))
  
}




## Deriving markers without logfc constraint, and allowing one other cell type to not be sig


# Main loop
all_markers <- list()
#initialise empty df 
all_markers_df <- data.frame()
celltypes <- unique(meta_data$assignments)
for (celltype in celltypes) {
  other_celltypes <- setdiff(celltypes, celltype)
  significant_genes <- compare_celltype(celltype, other_celltypes, design, fit, 0)
  print(paste("Cell type:", celltype, "Significant genes found:", nrow(significant_genes)))
  # Find genes significant in all but at most one comparison
  min_required <- length(other_celltypes)
  markers <- significant_genes %>% group_by(gene) %>% dplyr::filter(n() >= min_required-3) 
  print(paste("Markers after filtering for consistency in", celltype, ":", nrow(markers)))
  #only keep if log2fc is ALWAYS > 2
  markers <- markers %>% group_by(gene) %>% dplyr::filter(all(logFC > 0)) %>% dplyr::summarise(log2fc = mean(logFC), pval = exp(mean(log(adj.P.Val))), avg_expr = mean(AveExpr))
  #order by pval
  print(paste("Markers after filtering for log2fc > 2 in", celltype, ":", nrow(markers)))
  markers <- markers %>% dplyr::arrange(pval)
  #add col about cell types
  markers$celltype <- celltype
  
  
  #add markers to all_markers_df
  all_markers_df <- rbind(all_markers_df, markers)
  
  #remove dups
  markers <- unique(markers%>% pull(gene))
  all_markers[[celltype]] <- markers
}


all_markers_df <- all_markers_df %>% dplyr::arrange(pval)


write.csv(all_markers_df, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/relaxed_markers_mouse.csv")
