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
#remove where its less than 1.78 - log10 60

expression_table <- expression_table[,meta_data$Age_numeric >= log10(60)]
meta_data <- meta_data[meta_data$Age_numeric >= log10(60),]


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

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)

#make boxplot of colsums log y
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")


exp_table = expression_table[,meta_data$exptype == "sc"]
meta_data = meta_data[meta_data$exptype == "sc",]

dge = DGEList(counts = exp_table, group = meta_data$assignments)

# Create the design matrix
design <- model.matrix(~0 + assignments +  assignments:Age_numeric , data = meta_data)


#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 500,min.prop= 0.2)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)

#norm
dge <- calcNormFactors(dge, method = "TMM")

colnames(design)<-make.names(colnames(design))


param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~0 + assignments +  assignments:Age_numeric +  (1 | Author)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, meta_data, BPPARAM = param)

# fit dream model with contrasts
fit <- dream(vobjDream, form, meta_data)

head(coef(fit))

fit <- eBayes(fit)

top_genes <- topTable(fit, coef = "assignmentsStem_cells:Age_numeric" ,number = Inf, adjust.method="BH", lfc=0)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))


#how many sig
print(nrow(top_genes %>% dplyr::filter(adj.P.Val < 0.05 )))



background_genes <- rownames(dge)
#save to /Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/enrichment/aging_enrichment/aging_analysis_summary.csv
write.csv(background_genes, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/enrichment/aging_enrichment/sc_aging_background_genes_2_month.csv", row.names = FALSE)


#####
#Age effect
#####


cpdb <- read_csv("/Users/k23030440/epitome_code/epitome/data/gene_group_annotation/v_0.01/cpdb.csv")
#keep where category is TF and keep gene
tf_list <- cpdb %>% dplyr::filter(category == "TF") %>% dplyr::select(gene) %>% distinct()
tf_list <- tf_list$gene


print("Up TFs with time")
#rank it by adj.P.val
print(top_genes %>% dplyr::filter(genes %in% tf_list & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))
print("Down TFs with time")
print(top_genes %>% dplyr::filter(genes %in% tf_list & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))


print("Up with time")
#rank it by adj.P.val
print(top_genes %>% dplyr::filter( logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))
print("Down  with time")
print(top_genes %>% dplyr::filter( logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val))



#save these tfs up and down
tfs_up = top_genes %>% dplyr::filter(genes %in% tf_list & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
tfs_down = top_genes %>% dplyr::filter(genes %in% tf_list & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)

tfs_up

tfs_down

#save to "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
write.csv(tfs_up, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_tfs_up_2months.csv")
write.csv(tfs_down, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_tfs_down_2months.csv")

cpdb <- read_csv("/Users/k23030440/epitome_code/epitome/data/gene_group_annotation/v_0.01/cpdb.csv")
#keep where category is TF and keep gene
ligands  <- cpdb %>% dplyr::filter(category == "ligand") %>% dplyr::select(gene) %>% distinct()
ligands <- ligands$gene
receptors  <- cpdb %>% dplyr::filter(category == "receptor") %>% dplyr::select(gene) %>% distinct()
receptors <- receptors$gene

#ligands /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/CellChatDBmouse.csv
ligands_up = top_genes %>% dplyr::filter(genes %in% ligands & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
ligands_down = top_genes %>% dplyr::filter(genes %in% ligands & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
ligands_up
ligands_down


#receptors
receptors_up = top_genes %>% dplyr::filter(genes %in% receptors & logFC > 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
receptors_down = top_genes %>% dplyr::filter(genes %in% receptors & logFC < 0 & adj.P.Val < 0.05) %>% arrange(adj.P.Val)
receptors_up
receptors_down

write.csv(ligands_up, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_ligands_up_2months.csv")
write.csv(ligands_down, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_ligands_down_2months.csv")
write.csv(receptors_up, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_receptors_up_2months.csv")
write.csv(receptors_down, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sc_receptors_down_2months.csv")


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



coefs = c("assignmentsCorticotrophs:Age_numeric",
          "assignmentsEndothelial_cells:Age_numeric", "assignmentsGonadotrophs:Age_numeric",     
          "assignmentsImmune_cells:Age_numeric",      "assignmentsLactotrophs:Age_numeric",      
          "assignmentsMelanotrophs:Age_numeric",      "assignmentsMesenchymal_cells:Age_numeric",
          "assignmentsPituicytes:Age_numeric",        "assignmentsSomatotrophs:Age_numeric",     
          "assignmentsStem_cells:Age_numeric",        "assignmentsThyrotrophs:Age_numeric" )


# Create empty list to store results
all_results <- list()

# Extract cell types from coefficient names
get_cell_type <- function(coef_name) {
  # Remove "assignments" prefix and ".Age_numeric" suffix
  gsub("assignments(.*)\\.Age_numeric", "\\1", coef_name)
}

# Iterate through each coefficient
for(coef in coefs) {
  print(coef)
 
  top_genes <- topTable(fit, coef=coef,number = Inf,adjust.method="bonf", lfc=0.5)
  top_genes$genes <- rownames(top_genes)
  #print cell type
  print(get_cell_type(coef))
  #print number of significant genes
  print(nrow(top_genes %>% dplyr::filter(adj.P.Val < 0.05 )))
  
  # Add cell type column
  top_genes$cell_type <- get_cell_type(coef)
  
  # Store results
  all_results[[coef]] <- top_genes
}

# Combine all results into single dataframe
final_results <- do.call(rbind, all_results)

# Reset row names
rownames(final_results) <- NULL

#add col sig, and sign
final_results$sig <- ifelse(final_results$adj.P.Val < 0.05, "yes", "no")
final_results$sign <- ifelse(final_results$logFC > 0, "up", "down")



#save to /Users/k23030440/epitome_code/epitome/aging/v_0.02
version = "v_0.01"
write.csv(final_results, paste0("/Users/k23030440/epitome_code/epitome/data/aging/v_0.01/aging_genes_2_months.csv"))


final_results <- read.csv(paste0("/Users/k23030440/epitome_code/epitome/data/aging/v_0.01/aging_genes_2_months.csv"), header = TRUE, stringsAsFactors = FALSE)
final_results_all_ages <- read.csv(paste0("/Users/k23030440/epitome_code/epitome/data/aging/v_0.01/aging_genes.csv"), header = TRUE, stringsAsFactors = FALSE)
final_results$sig <- ifelse(final_results$adj.P.Val < 0.05, "yes", "no")
final_results_all_ages$sig <- ifelse(final_results_all_ages$adj.P.Val < 0.05, "yes", "no")
final_results$direction <- ifelse(final_results$logFC > 0, "up", "down")
final_results_all_ages$direction <- ifelse(final_results_all_ages$logFC > 0, "up", "down")
table(final_results$cell_type, final_results$sig, final_results$direction)
table(final_results_all_ages$cell_type, final_results_all_ages$sig, final_results_all_ages$direction)




library(dplyr)
library(stringr)
library(ggplot2)

# keep significant genes only
sig_2mo <- final_results %>%
  filter(sig == "yes") %>%
  mutate(
    cell_type = str_replace(cell_type, "assignments", ""),
    cell_type = str_replace(cell_type, ":Age_numeric", ""),
    analysis = "2mo"
  )

sig_all <- final_results_all_ages %>%
  filter(sig == "yes") %>%
  mutate(
    analysis = "all_ages"
  )

sig_combined <- full_join(
  sig_2mo,
  sig_all,
  by = c("genes", "cell_type", "direction"),
  suffix = c("_2mo", "_all")
) %>%
  mutate(
    category = case_when(
      !is.na(adj.P.Val_2mo) & !is.na(adj.P.Val_all) ~ "Both",
      !is.na(adj.P.Val_2mo) &  is.na(adj.P.Val_all) ~ "2mo_only",
      is.na(adj.P.Val_2mo) & !is.na(adj.P.Val_all) ~ "All_ages_only"
    )
  )

plot_df <- sig_combined %>%
  group_by(cell_type, direction, category) %>%
  summarise(n_genes = n(), .groups = "drop")

plot_df$category <- factor(
  plot_df$category,
  levels = c("2mo_only", "Both", "All_ages_only")
)


plot_stacked <- function(df, dir_label, title_suffix) {
  ggplot(
    df %>% filter(direction == dir_label),
    aes(x = cell_type, y = n_genes, fill = category)
  ) +
    geom_bar(stat = "identity") +
    labs(
      x = "Cell type",
      y = "Number of DE genes",
      fill = "Analysis source",
      title = paste(title_suffix, dir_label, "regulated genes")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

# Upregulated
plot_up <- plot_stacked(
  plot_df,
  dir_label = "up",
  title_suffix = "Differentially expressed"
)

# Downregulated
plot_down <- plot_stacked(
  plot_df,
  dir_label = "down",
  title_suffix = "Differentially expressed"
)

plot_up
plot_down


#get examples from Stem_cells where its only in 2 months
sig_combined %>%
  filter(cell_type == "Stem_cells", direction == "up", category == "2mo_only") %>%
  arrange(adj.P.Val_2mo)


sig_combined %>%
  filter(cell_type == "Stem_cells", direction == "up", category == "All_ages_only") %>%
  arrange(adj.P.Val_2mo)


sig_combined %>%
  filter(cell_type == "Stem_cells", direction == "up", category == "Both") %>%
  arrange(adj.P.Val_2mo)

