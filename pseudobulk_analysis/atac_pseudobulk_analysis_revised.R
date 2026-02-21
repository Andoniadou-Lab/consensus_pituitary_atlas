#loading libraries
library(DESeq2)
library(limma)
library(edgeR)
library(tidyverse)
library(org.Mm.eg.db)
library(reticulate)
library(Seurat)


# Specify the path to the .h5ad file
#h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pb_h5ad.h5ad"
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pb_h5ad_0904.h5ad"
# Read the .h5ad file using reticulate
anndata <- reticulate::import("anndata")
py_index <- reticulate::import("pandas")$Index
adata <- anndata$read_h5ad(h5ad_file)

multiome_path = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pdata_assigned_0828.h5ad"

multiome <- t(anndata$read_h5ad(multiome_path)$X)

multiome_metadata = as.data.frame(anndata$read_h5ad(multiome_path)$obs)

#cell_type_final rename to cell_type, GEO to sample
colnames(multiome_metadata)[colnames(multiome_metadata) == "cell_type_final"] <- "cell_type"
colnames(multiome_metadata)[colnames(multiome_metadata) == "GEO"] <- "sample"

vars_multi = py_to_r(anndata$read_h5ad(multiome_path)$var_names$to_list())

# Convert the data to Seurat object
# Seurat expects the count matrix to be in column-major order (genes x cells)
counts <- t(adata$X) # Transpose to convert to column-major order
meta_data <- as.data.frame(adata$obs) # Convert obs to a data frame
colnames(meta_data)[colnames(meta_data) == "GEO"] <- "sample"
vars <- adata$var
vars <-  py_to_r(adata$var_names$to_list())
vars

#whats the overlap of vars and vars_multi
print(length(intersect(vars, vars_multi)))

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = counts, meta.data = meta_data)
rownames(seurat_object) <- vars
seurat_object_multi<- CreateSeuratObject(counts = multiome, meta.data = multiome_metadata)
rownames(seurat_object_multi) <- vars_multi

#merge all
all_seurat <- merge(seurat_object, y = c(seurat_object_multi),
                    add.cell.ids = c("all","multi"),
                    project = "atac_pseudobulk")

all_seurat<- JoinLayers(all_seurat)

#print total counts for all cells separately
print(colSums(all_seurat@assays$RNA$counts))
      
seurat_object<- all_seurat
rownames(seurat_object)

#save as /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/consensus_chromatin_landscape.csv
write.csv(vars, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/consensus_chromatin_landscape.csv", row.names = FALSE)
write.csv(vars, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/additional_downloads/consensus_chromatin_landscape/v_0.01/consensus_chromatin_landscape.csv", row.names = FALSE)
#normalize
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1000000)

root = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/"
figs_folder = paste(root,"Figures/",sep="")


meta_data = seurat_object@meta.data
meta_data$index <- 1:nrow(meta_data)
#rename sample to GEO
meta_data$GEO = meta_data$sample
print(length(unique(meta_data$Author)))
print(length(unique(meta_data$SRA_ID)))


meta2 = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/curation/v_0.01/cpa.csv")
#rename GEO to sample
colnames(meta2)[colnames(meta2) == "GEO"] <- "sample"
#add missing columns to meta_data from meta2. shared columns is GEO. keep the same row dimensions
#iterate through meta_data
# First make meta2 have only first occurrence of each GEO
meta2_first <- meta2[!duplicated(meta2$sample), ]

#metadata only keep sample column
meta_data = meta_data[,c("sample","index")]
# Then merge will maintain row dimensions
meta_data <- merge(meta_data, meta2_first, by="sample", all.x=TRUE)
#reorder based on index
meta_data <- meta_data[order(meta_data$index),]

seurat_object <- AddMetaData(seurat_object, meta_data)

#print unique SRA_ID
print(length(unique(meta_data$SRA_ID)))

meta_data %>%
  distinct(SRA_ID, Comp_sex) %>%  
  count(Comp_sex, name = "num_SRA_IDs")


meta_data = seurat_object@meta.data

expression_table = seurat_object@assays$RNA$counts
#make a hist of colsums
hist(colSums(expression_table), breaks = 100)
#print median value of the table :200, :20000
print(median(expression_table[1:200,1:200]))

#print hist of cells meta_data$psbulk_n_cells
hist(as.numeric(meta_data$psbulk_n_cells), breaks = 100)

##########################

library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$cell_type)
meta_data$exptype <- make.names(meta_data$Modality)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)


#remove where assignment is Blood
expression_table <- expression_table[,meta_data$assignments != "Erythrocytes"]
meta_data <- meta_data[meta_data$assignments != "Erythrocytes",]

#where Modality is multi_rna turn to multi_atac
meta_data$exptype[meta_data$Modality == "multi_rna"] <- "multi_atac"



# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$Modality)
meta_data$Sex <- make.names(meta_data$Comp)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)
meta_data$Sex <- gsub("[^A-Za-z0-9_]", "", meta_data$Sex)

print(meta_data$exptype)
#add a new column called index
meta_data$index <- 1:nrow(meta_data)
#make index the rownames
rownames(meta_data) <- meta_data$index

#divide values in each column by respective value in meta_data meta_data$psbulk_n_cells
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")

# Convert to numeric
meta_data$psbulk_n_cells <- as.numeric(as.character(meta_data$psbulk_n_cells))



expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)


#make boxplot of colsums log y
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")

#where multi_rna changer to multi_atac
meta_data$exptype[meta_data$exptype == "multi_rna"] <- "multi_atac"


dge = DGEList(counts = expression_table, group = meta_data$assignments)
dge <- calcNormFactors(dge, method = "TMMwsp")


# Create the design matrix
design <- model.matrix(~ 0 + assignments , data = meta_data) 

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 300, min.prop= 0.1)
print(sum(keep_genes))
#pull out rownames starting with mt or rpl or rps
dge <- dge[keep_genes, ]
dim(dge)

fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

###
# Make column names valid
colnames(design) <- make.names(colnames(design), unique = TRUE)

print(head(coef(fit)))
#turn coef(fit) into a dataframe
coef_df = as.data.frame(coef(fit))

#save 
write.csv(as.data.frame(coef(fit)), "/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/coef.csv")

coef_df = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/coef.csv")

dim(coef_df)
coef_df$gene = rownames(coef_df)

# Create contrast
contrast <- makeContrasts(
  #assignmentsCorticotrophs - assignmentsMelanotrophs, 
  assignmentsCorticotrophs - assignmentsMelanotrophs,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, trend = TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf, sort.by = "P", lfc = 0.5)
top_genes$genes <- rownames(top_genes)
print(head(top_genes))

#print how many are adj.P.val < 0.05
print(sum(top_genes$adj.P.Val < 0.05))

#make scatterplot of log2fc and pval
plot(top_genes$logFC, -log10(top_genes$P.Value), alpha=0.1,xlab = "log2 fold change", ylab = "-log10 p-value", main = "Volcano plot of differential expression", pch = 20, col = ifelse(top_genes$adj.P.Val < 0.05, "red", "black"))


#print total number of peaks
print(nrow(top_genes))

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
    
    sig_genes <- topTable(fit2, number = Inf, lfc = 0.5) %>%
      rownames_to_column("gene") %>%
      dplyr::filter(adj.P.Val < 0.05 & logFC > log2fc) %>%
      mutate(group1 = celltype,
             group2 = other)
    
    significant_genes <- bind_rows(significant_genes, sig_genes)
  }
  
  return(significant_genes)
}

# Perform comparisons for all cell types
all_comparisons <- data.frame()

library(tidyverse)
for (celltype in celltypes) {
  other_celltypes <- celltypes[celltypes != celltype]
  comparison_results <- compare_celltype(celltype, other_celltypes, design, fit)
  all_comparisons <- bind_rows(all_comparisons, comparison_results)
}

# Display the first few rows of the results
print(head(all_comparisons))


#I have a table called all_comparisons with columns "gene"      "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"        "group1"    "group2"
#Below I defined different groupings of comparisons. Here I am trying to find the genes that are significant in all comparisons of a grouping, and in the same direction.
#Give me the corresponding gene_lists, and whether they are up or down regulated in each comparison

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

grouping_lists = list(grouping_1, grouping_2, grouping_3, grouping_4, grouping_5, grouping_6, grouping_7, grouping_8)

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
               stringsAsFactors = FALSE)
  }
  
  # Get significant genes for all comparisons and stitch to a single df
  sig_genes_list <- data.frame(gene = character(), direction = character(), stringsAsFactors = FALSE)
  for (i in 1:length(comps_left)) {
    new_genes <- get_sig_genes(data.frame(group1 = comps_left[i], group2 = comps_right[i]))
    sig_genes_list <- rbind(sig_genes_list, new_genes)
  }
  
  library(dplyr)
  # Only keep genes that occur in all comparisons - meaning n_comps times
  sig_genes_list <- sig_genes_list %>% 
    group_by(gene) %>%
    dplyr::filter(n() == n_comps) %>%
    ungroup()
  
  # Only keep genes where the direction is the same in all comparisons
  sig_genes_list <- sig_genes_list %>% 
    group_by(gene) %>%
    dplyr::filter(n_distinct(direction) == 1) %>%
    ungroup()
  
  #add column mean_log2fc
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% mutate(mean_log2fc = mean(log2fc))
  #add geom_mean_adj_pval
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% mutate(geom_mean_adj_pval = exp(mean(log(pvalue))))
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

write_csv(all_results, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/grouping_lineage_markers.csv")
write_csv(all_results, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/heatmap/v_0.01/atac_grouping_lineage_markers.csv")
write.csv(all_results, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/markers/v_0.01/atac_grouping_lineage_markers.csv")


all_results <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/grouping_lineage_markers.csv")


#keep obly sig
all_results <- all_results[all_results$geom_mean_adj_pval<0.05,]
#how many per grouping
table(all_results$grouping, all_results$direction)

for (i in 1:8) {
  grouping_results <- all_results %>% dplyr::filter(grouping == paste0("grouping_", i))
  cat("\nGrouping", i, "results:\n")
  cat("Number of significant peaks:", nrow(grouping_results), "\n")
}


library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)

grouping1_up = all_results[all_results$grouping=="grouping_1",]
#also direction up
grouping1_up = grouping1_up[grouping1_up$direction=="down",]
#make mean log2fc positive > 2
#grouping1_up = grouping1_up[grouping1_up$mean_log2fc>2,]
rest = all_results[all_results$grouping!="grouping_1",]


##############
#Signac enrichment
##############

library(Signac)
atac_assay <- CreateChromatinAssay(seurat_object[["RNA"]]$counts,sep = c(":", "-"),)

atac_object <- seurat_object
atac_object[["ATAC"]]<- atac_assay
DefaultAssay(atac_object) <- "ATAC"
atac_object[["RNA"]] <- NULL

library(JASPAR2022)

pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

seqlevelsStyle(BSgenome.Mmusculus.UCSC.mm10) <- "UCSC"

# add motif information
atac_object <- AddMotifs(
  object = atac_object,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm,
  p.cutoff = 1e-05, w = 10
)



target_peaks <- all_results[all_results$grouping == "grouping_2" & 
                              all_results$direction == "up", ]

# Skip if no peaks found
if(nrow(target_peaks) == 0) {
  message(sprintf("No peaks found for %s (%s)", "grouping_1", "down"))
  return(NULL)
}

#grouping1_up$gene change : to -
target_peaks$gene <- gsub(":", "-", target_peaks$gene)

all_target_peaks <- all_results#[all_results$grouping == "grouping_2", ]
all_target_peaks$gene <- gsub(":", "-", all_target_peaks$gene)

#keep only target_peaks
atac_object_subset = atac_object[rownames(atac_object) %in% all_target_peaks$gene]

enriched.motifs <- FindMotifs(
  object = atac_object_subset,
  features = target_peaks$gene
)





##########################################
##########################################
##################### Quick analysis of those peaks with JUN sites (MA1144.1, MA1134.1, MA1137.1, MA1130.1)
#create peaks_df with cols peak, mean, var, JUN site
# Extract peak names
peak <- rownames(atac_object)

# Calculate mean and variance for each peak row in ATAC data at once
means <- rowMeans(atac_object[["ATAC"]]$data)
variances <- apply(atac_object[["ATAC"]]$data, 1, var)

# Filter JUN motif columns and determine if each row has any JUN site
motif_df_1 <- atac_object@assays$ATAC@motifs@data[, c("MA1144.1", "MA1134.1", "MA1137.1", "MA1130.1")]
JUN_sites <- rowSums(motif_df_1) > 0

# Create the data frame in one go
peaks_df <- data.frame(
  peak = peak,
  mean = means,
  var = variances,
  JUN_site = JUN_sites,
  stringsAsFactors = FALSE
)

#make scatter of mean and var and color by JUN site. log x and y
ggplot(peaks_df, aes(x = log(mean), y = log(var), color = JUN_site)) +
  geom_point(alpha = 0.05) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(x = "log(mean)", y = "log(var)")

#plot only red
ggplot(peaks_df[peaks_df$JUN_site, ], aes(x = log(mean), y = log(var))) +
  geom_point(color = "red",alpha = 0.05) +
  theme_minimal() +
  labs(x = "log(mean)", y = "log(var)")

#plot only black
ggplot(peaks_df[!peaks_df$JUN_site, ], aes(x = log(mean), y = log(var))) +
  geom_point(color = "black",alpha = 0.05) +
  theme_minimal() +
  labs(x = "log(mean)", y = "log(var)")
  


##########################################
##########################################
##########################################


# Function to perform motif analysis for a specific grouping and direction
perform_motif_analysis <- function(all_results, grouping_name, direction, output_dir) {
  # Filter for specific grouping and direction
  
  target_peaks <- all_results[all_results$grouping == grouping_name & 
                                all_results$direction == direction, ]
  
  # Skip if no peaks found
  if(nrow(target_peaks) == 0) {
    message(sprintf("No peaks found for %s (%s)", grouping_name, direction))
    return(NULL)
  }
  #grouping1_up$gene change : to -
  target_peaks$gene <- gsub(":", "-", target_peaks$gene)
  
  
  
  all_target_peaks <- all_results[all_results$grouping == grouping_name, ]
  all_target_peaks$gene <- gsub(":", "-", all_target_peaks$gene)
  #keep only target_peaks
  atac_object_subset = atac_object[rownames(atac_object) %in% all_target_peaks$gene]
  
  #atac_object_ only keep peaks that start with chr
  atac_object_ = atac_object[grep("^chr", rownames(atac_object)), ]
  target_peaks$gene = target_peaks$gene[grep("^chr", target_peaks$gene)]
  
  enriched.motifs <- FindMotifs(
    object = atac_object_,
    features = target_peaks$gene,
    background = 100000
  )
  enriched.motifs
  output_file = file.path(output_dir, sprintf("enrichment_results_%s_%s.rds", grouping_name, direction))
    # Save results
  saveRDS(enriched.motifs, output_file)
  write.csv(
    enriched.motifs,
        file.path(output_dir, 
                  sprintf("enrichment_results_%s_%s.csv", 
                          grouping_name, 
                          direction)),
        row.names = FALSE
      )
  
  
  write.csv(
    enriched.motifs,
    file.path("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/v_0.01/", 
              sprintf("enrichment_results_%s_%s.csv", 
                      grouping_name, 
                      direction)),
    row.names = FALSE
  )
  
    
  return(enriched.motifs)
}

# Create output directory
output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_analysis_results"
dir.create(output_dir, showWarnings = FALSE)

# Get unique groupings
groupings <- unique(all_results$grouping)
directions <- c("up", "down")

# Initialize list to store all results
all_enrichment_results <- list()

# Loop through each grouping and direction
for(grouping in groupings) {
  for(direction in directions) {
    message(sprintf("Processing %s (%s)...", grouping, direction))
    
    result <- perform_motif_analysis(all_results, 
                                     grouping, 
                                     direction, 
                                     output_dir)
    
    if(!is.null(result)) {
      all_enrichment_results[[paste(grouping, direction, sep="_")]] <- result
    }
  }
}

# Save complete results object
saveRDS(all_enrichment_results, file.path(output_dir, "all_motif_results.rds"))

# Save the motif matrix
matrix = atac_object@assays$ATAC@motifs@data
write.table(rownames(matrix),
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_matrix_rows.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(rownames(matrix),
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/heatmap/v_0.01/motif_matrix_rows.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


write.table(colnames(matrix),
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_matrix_cols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(matrix),
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/heatmap/v_0.01/motif_matrix_cols.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

library(Matrix)
writeMM(matrix, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_matrix.mtx")
writeMM(matrix, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/heatmap/v_0.01/motif_matrix.mtx")
#Saving normalised data



# Preprocessing script to normalize and save data
library(edgeR)

# Apply normalization
normalize_expression <- function(expression_table_original, dge_original) {
  norm_factors <- dge_original$samples$norm.factors
  lib_sizes <- dge_original$samples$lib.size
  
  eff_lib_sizes <- norm_factors * lib_sizes
  
  normalized_expression <- sweep(expression_table_original, 2, eff_lib_sizes, FUN = "/")
  
  # Multiply back by 10^6
  normalized_expression <- t(t(normalized_expression) * 10^6)
  
  return(normalized_expression)
}

# Normalize the data
normalized_expression <- normalize_expression(expression_table, dge)

# Save the normalized data and metadata
save(normalized_expression, meta_data, 
     file = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/atac_normalized_data.RData")


#save meta_data as well
write.csv(meta_data, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/v_0.01/atac_meta_data.csv")


write.table(rownames(normalized_expression),
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/v_0.01/accessibility_features.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(colnames(normalized_expression),
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/v_0.01/accessibility_columns.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


library(Matrix)
normalized_expression <- as(normalized_expression, "CsparseMatrix")
writeMM(normalized_expression, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/v_0.01/normalized_data.mtx")


# Create summary table
summary_table <- data.frame()

for(name in names(all_enrichment_results)) {
  print(name)
  result <- all_enrichment_results[[name]]
    # Get top 10 motifs
    top_motifs <- head(result, 10)
    top_motifs$analysis <- name
    summary_table <- rbind(summary_table, top_motifs)
  
}

#also make a csv for all results
all_results_df = data.frame()
for(name in names(all_enrichment_results)) {
  print(name)
  result <- all_enrichment_results[[name]]
  result$analysis <- name
  all_results_df <- rbind(all_results_df, result)
}

# Save summary table
write.csv(summary_table, 
          file.path(output_dir, "motif_analysis_summary.csv"), 
          row.names = FALSE)

write.csv(all_results_df,
          file.path(output_dir, "all_motif_results.csv"),
          row.names = FALSE)


####################
####################
####################
####################
####################
####################
####################
####################
####################
####################
#Chromvar
rm(seurat_object)
rm(expression_table)
rm(atac_object_subset)
rm(atac_assay)



#save atac_object rds
saveRDS(atac_object, 
        "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/atac_object_0923.rds")

atac_object <- readRDS(
  "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/atac_object_0923.rds"
)

meta_data = atac_object@meta.data
print(colSums(atac_object@assays$ATAC$counts))

#iterate through unique SRA ids in meta_data. For each, make as subset and run chromvar. Then concat the resulting seurat objects back together in the end.
#print value counts for SRA_IDs
print(table(meta_data$SRA_ID))

atac_object<- RunChromVAR(
object = atac_object,
genome = BSgenome.Mmusculus.UCSC.mm10)



motifs = rownames(atac_object[["chromvar"]])
new_rownames = c()
for (motif in motifs) {
  #in atac_object@assays$ATAC@motifs@motif.names find the entry correspondign under motif
  motif_name = atac_object@assays$ATAC@motifs@motif.names[motif]
  new_rownames = c(new_rownames, paste0(motif, "_", motif_name))
}

matrix_vals = as.matrix(atac_object[["chromvar"]]@data)
matrix_vals = atac_object[["chromvar"]]@data

write.table(new_rownames, 
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/chromvar/v_0.01/chromvar_features.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(colnames(matrix_vals),
            "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/chromvar/v_0.01/chromvar_columns.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

library(Matrix)
matrix_vals <- as(as.matrix(atac_object[["chromvar"]]@data), "dgCMatrix")
writeMM(matrix_vals, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/chromvar/v_0.01/normalized_data.mtx")


#save meta_data as well
write_csv(atac_object@meta.data, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/chromvar/v_0.01/meta_data.csv")



#################
#################
#extracting motif locations for epitome
#################

#motif_names = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/chromvar/v_0.01/chromvar_features.txt"
motif_names = read.table("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/chromvar/v_0.01/chromvar_features.txt", header = FALSE, stringsAsFactors = FALSE)
motif_names = motif_names$V1
#choose values before _
motifs = sub("_.*", "", motif_names)
motifs
#add col
for (i in 1:length(atac_object[["ATAC"]]@motifs@positions)) {
  if (i %% 100 == 0) {
    print(i)
  }
  if (i == 1) {
    final_df = as.data.frame(atac_object[["ATAC"]]@motifs@positions[i])
    #turn to df
    final_df = as.data.frame(final_df)
    motif_name = unname(atac_object@assays$ATAC@motifs@motif.names[i])[[1]]
    motif = motifs[i]
    new_motif_name = paste0(motif, "_", motif_name)
    #add to col "motif"
    final_df$motif = new_motif_name
  }
  
  df = as.data.frame(atac_object[["ATAC"]]@motifs@positions[i])
  #turn to df
  df = as.data.frame(df)
  motif_name = unname(atac_object@assays$ATAC@motifs@motif.names[i])[[1]]
  motif = motifs[i]
  new_motif_name = paste0(motif, "_", motif_name)
  #add to col "motif"
  df$motif = new_motif_name
  
  final_df = rbind(final_df, df)
}

final_df <- final_df[!duplicated(final_df[, c("start","seqnames", "motif")]), ]

#save final_df

write.csv(final_df, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/v_0.01/atac_motif_data.csv")

#remove redundant combinations of start and motif
#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/accessibility/atac_meta_data.csv

# Load required libraries
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicFeatures)

# Load TxDb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

### Process Exons ###
# Extract exons grouped by gene
exons_by_gene <- exonsBy(txdb, by = "gene")
# Unlist to create a flat GRanges object for exons
exons <- unlist(exons_by_gene)
# Add gene IDs as metadata
exons$gene_id <- names(exons)
# Map gene IDs to gene names
exon_gene_ids <- mapIds(org.Mm.eg.db, 
                        keys = exons$gene_id, 
                        column = "SYMBOL", 
                        keytype = "ENTREZID", 
                        multiVals = "first")
mcols(exons)$gene_name <- exon_gene_ids

# Convert to a data frame without using gene IDs as row names
exon_df <- as.data.frame(exons, row.names = NULL)
exon_df$feature_type <- "exon"

### Process Introns ###
# Extract introns grouped by transcript
introns_by_gene <- intronsByTranscript(txdb)
# Unlist to create a flat GRanges object for introns
introns <- unlist(introns_by_gene)
# Map transcript IDs to gene names
intron_gene_ids <- mapIds(org.Mm.eg.db, 
                          keys = names(introns), 
                          column = "SYMBOL", 
                          keytype = "ENTREZID", 
                          multiVals = "first")
mcols(introns)$gene_name <- intron_gene_ids

# Convert to a data frame without using transcript IDs as row names
intron_df <- as.data.frame(introns, row.names = NULL)
intron_df$feature_type <- "intron"

### Combine Exons and Introns ###
# Add unique identifiers if needed
exon_df$feature_id <- paste0("exon_", seq_len(nrow(exon_df)))
intron_df$feature_id <- paste0("intron_", seq_len(nrow(intron_df)))

# Combine exon and intron data frames
#keep thsese from exon_df seqnames   start     end width strand gene_name feature_type feature_id
exon_df <- subset(exon_df, select = c("seqnames", "start", "end", "width", "strand", "gene_name", "feature_type", "feature_id"))
annotation_df <- rbind(exon_df)

#save annotation_df to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/accessibility/annotation.csv
write.csv(annotation_df, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/accessibility/v_0.01/annotation.csv", row.names = FALSE)



#find Sox2 in annotation_df
sox2 <- subset(annotation_df, gene_name == "Sox2")

sox2




##################
##################
##################
#Looking at peaks changing with sex

# Specify the path to the .h5ad file
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pb_h5ad_0904.h5ad"

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


#normalize
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1000000)
VlnPlot(seurat_object, features = c("chr1:89456336-89456837"),group.by = "cell_type", pt.size = 0.1)

root = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_one/Proposal/pituitary atlas/"
figs_folder = paste(root,"Figures/",sep="")


meta_data = seurat_object@meta.data
meta_data$index <- 1:nrow(meta_data)
#rename sample to GEO
meta_data$GEO = meta_data$sample
print(length(unique(meta_data$Author)))
print(length(unique(meta_data$SRA_ID)))



meta2 = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/epitome/data/curation/v_0.01/cpa.csv")
#add missing columns to meta_data from meta2. shared columns is GEO. keep the same row dimensions
#iterate through meta_data
# First make meta2 have only first occurrence of each GEO
meta2_first <- meta2[!duplicated(meta2$GEO), ]

# Then merge will maintain row dimensions
meta_data <- merge(meta_data, meta2_first, by="GEO", all.x=TRUE)
#reorder based on index
meta_data <- meta_data[order(meta_data$index),]

seurat_object <- AddMetaData(seurat_object, meta_data)
#remove obs that are Normal != 1
seurat_object 
print(dim(seurat_object))


meta_data = seurat_object@meta.data

expression_table = seurat_object@assays$RNA$counts

##########################

library(limma)
library(edgeR)
# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$cell_type)
meta_data$exptype <- make.names(meta_data$Modality)


# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)


#remove where assignment is Blood
expression_table <- expression_table[,meta_data$assignments != "Erythrocytes"]
meta_data <- meta_data[meta_data$assignments != "Erythrocytes",]


#where Modality is multi_rna turn to multi_atac
meta_data$exptype[meta_data$Modality == "multi_rna"] <- "multi_atac"

meta_data$Age_numeric <- as.numeric(as.character(meta_data$Age_numeric))
meta_data$Comp<- as.numeric(as.character(meta_data$Comp_sex))





# Ensure valid column names for assignments and exptype
meta_data$assignments <- make.names(meta_data$assignments)
meta_data$exptype <- make.names(meta_data$Modality)
meta_data$Sex <- make.names(meta_data$Comp)

# Change 'multi_rna' to 'sn' in exptype
#meta_data$exptype[meta_data$exptype == "multi_rna"] <- "sn"

# Remove any remaining problematic characters
meta_data$assignments <- gsub("[^A-Za-z0-9_]", "", meta_data$assignments)
meta_data$exptype <- gsub("[^A-Za-z0-9_]", "", meta_data$exptype)
meta_data$Sex <- gsub("[^A-Za-z0-9_]", "", meta_data$Sex)

#print all values
print(meta_data$exptype)
#add a new column called index
meta_data$index <- 1:nrow(meta_data)
#make index the rownames
rownames(meta_data) <- meta_data$index

#divide values in each column by respective value in meta_data meta_data$psbulk_n_cells
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")

# Convert to numeric
meta_data$psbulk_n_cells <- as.numeric(as.character(meta_data$psbulk_n_cells))

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)

#make boxplot of colsums log y
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")

#where multi_rna changer to multi_atac
meta_data$exptype[meta_data$exptype == "multi_rna"] <- "multi_atac"


dge = DGEList(counts = expression_table, group = meta_data$assignments)
#dge$samples = rbind(dge_atac$samples,dge_multi$samples)
dge <- calcNormFactors(dge, method = "TMMwsp")

# Create the design matrix
design <- model.matrix(~ 0 + assignments + assignments:Sex , data = meta_data)

#make name correct with make.names(
colnames(design) <- make.names(colnames(design))
#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- rowSums(dge$counts > 1) >= 0.01 * ncol(dge$counts) & rowSums(dge$counts) >= 1000

#pull out rownames starting with mt or rpl or rps
dge <- dge[keep_genes, ]
dim(dge)


fit <- voom(dge, design, plot=TRUE)

fit <- lmFit(fit, design)

colnames(fit)



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
top_genes <- topTable(diff_expression, number = Inf,  adjust.method = "BH", lfc = 1)
top_genes$genes <- rownames(top_genes)
#how many sig
nrow(top_genes[top_genes$adj.P.Val < 0.05, ])

library(clipr)
# Write the genes to clipboard, one per line
up_genes <- top_genes %>% dplyr::filter(logFC > 0 & adj.P.Val < 0.05)
clipr::write_clip(paste(up_genes$gene, collapse = "\n"))
down_genes <- top_genes %>% dplyr::filter(logFC < 0 & adj.P.Val < 0.05)
clipr::write_clip(paste(down_genes$gene, collapse = "\n"))



# Initialize a list to store results
sex_specific_regions <- list()
cell_types <- unique(meta_data$assignments)

# Print column names to debug
print("Design matrix column names:")
print(colnames(design))

results_list <- list()

for (cell_type in cell_types) {
  tryCatch({
    # Create contrast name - be explicit about the interaction term
    contrast_name <- paste0("assignments", cell_type, ".SexX1")
    
    # Check if the contrast exists in the design matrix
    if(!contrast_name %in% colnames(design)) {
      cat(sprintf("Skipping %s: contrast not found in design matrix\n", cell_type))
      next
    }
    
    # Create contrast matrix
    contrast <- makeContrasts(
      sex_effect = contrast_name,
      levels = colnames(design)
    )
    
    # Compute differential expression
    diff_expr <- contrasts.fit(fit, contrast)
    diff_expr <- eBayes(diff_expr, robust=TRUE)
    
    # Get results
    results <- topTable(diff_expr, number=Inf, lfc = 0.5)
    #remove regions that start with chrY
    #results <- results[!grepl("chrY",rownames(results)),]
    
    # Store significant regions
    male_regions <- rownames(results[results$adj.P.Val < 0.05 & results$logFC > 0,])
    female_regions <- rownames(results[results$adj.P.Val < 0.05 & results$logFC < 0,])
    
    sex_specific_regions[[cell_type]] <- list(
      male = male_regions,
      female = female_regions
    )
    
    cat(sprintf("%s: %d male-specific, %d female-specific regions\n", 
                cell_type, length(male_regions), length(female_regions)))
    
    results$cell_type <- cell_type
    results_list[[cell_type]] <- results
    
  }, error = function(e) {
    cat(sprintf("Error processing %s: %s\n", cell_type, e$message))
  })
}

results_list <- do.call(rbind, results_list)
#split rowname at . to celltype and peak
results_list$cell_type <- sapply(strsplit(rownames(results_list), "\\."), `[`, 1)
results_list$peak <- sapply(strsplit(rownames(results_list), "\\."), `[`, 2)
#keep adj.P.Val < 0.05
results_list <- results_list[results_list$adj.P.Val < 0.05,]
lactotrophs_up <- results_list[results_list$logFC > 0 & results_list$cell_type == "Lactotrophs",]
lactotrophs_down <- results_list[results_list$logFC < 0 & results_list$cell_type == "Lactotrophs",]
gonadotrophs_up <- results_list[results_list$logFC > 0 & results_list$cell_type == "Gonadotrophs",]
gonadotrophs_down <- results_list[results_list$logFC < 0 & results_list$cell_type == "Gonadotrophs",]









# Create plotting data
plot_data <- data.frame(
  cell_type = names(sex_specific_regions),
  male_regions = sapply(sex_specific_regions, function(x) length(x$male)),
  female_regions = sapply(sex_specific_regions, function(x) length(x$female))
)

# Create long format data
plot_data_long <- rbind(
  data.frame(
    cell_type = plot_data$cell_type,
    count = -plot_data$female_regions,  # Negative for female to show on left
    sex = "Female",
    label = plot_data$female_regions    # Positive label
  ),
  data.frame(
    cell_type = plot_data$cell_type,
    count = plot_data$male_regions,
    sex = "Male",
    label = plot_data$male_regions
  )
)

# Remove empty labels
plot_data_long$label[plot_data_long$label == 0] <- ""

# Add spacing for male labels
plot_data_long$label <- ifelse(plot_data_long$sex == "Male", 
                               paste0("   ", plot_data_long$label),
                               plot_data_long$label)

#where 0 remove
plot_data_long <- plot_data_long[plot_data_long$count != 0, ]

# Calculate total differentially accessible regions per cell type and reorder
plot_data_long <- plot_data_long %>%
  group_by(cell_type) %>%
  mutate(total_diff_acc = sum(abs(count))) %>%
  ungroup() %>%
  arrange(desc(total_diff_acc)) %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

# Add fixed label positions
plot_data_long <- plot_data_long %>%
  mutate(label_position = ifelse(count > 0, 5000, -5000))

# Create the plot
plt = ggplot(plot_data_long, aes(x = cell_type, y = count, fill = sex)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = label, y = label_position),
            color = "black",
            size = 5) +
  scale_fill_manual(values = c("Female" = "#FFA500", "Male" = "#63B3ED")) +
  scale_y_continuous(
    labels = abs,
    limits = c(-max(abs(plot_data_long$count)) * 1.2, max(abs(plot_data_long$count)) * 1.2)
  ) +
  coord_flip() +
  labs(
    title = "Sex-biased Chromatin Accessibility Across Cell Types",
    x = "Cell Type",
    y = "Number of Sex-biased Regions",
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

# Save the plot
fold = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/"
ggsave(plt, filename = paste0(fold,"sex_biased_accessibility.png"), width = 6, height = 5.5, dpi = 300)
ggsave(plt, filename = paste0(fold,"sex_biased_accessibility.svg"), width = 6, height = 5.5, dpi = 300)

plt



##############
#Signac enrichment
##############
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)

library(Signac)
atac_assay <- CreateChromatinAssay(seurat_object[["RNA"]]$counts,sep = c(":", "-"),)

atac_object <- seurat_object
atac_object[["ATAC"]]<- atac_assay
DefaultAssay(atac_object) <- "ATAC"
atac_object[["RNA"]] <- NULL


library(JASPAR2022)

pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

seqlevelsStyle(BSgenome.Mmusculus.UCSC.mm10) <- "UCSC"

library(EnsDb.Mmusculus.v79)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# Change to UCSC style to match the genome and data
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# Add the gene information to the ATAC object
Annotation(atac_object) <- annotations




# add motif information
atac_object <- AddMotifs(
  object = atac_object,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm,
  p.cutoff = 1e-05, w = 10
)


# Function to perform motif analysis for a specific grouping and direction
perform_motif_analysis <- function(sex_specific_regions, grouping_name, direction, output_dir) {
  # Filter for specific grouping and direction
  
  target_peaks <- sex_specific_regions[[grouping_name]][[direction]]
  if(is.null(target_peaks)) {
    message(sprintf("No peaks found for %s (%s)", grouping_name, direction))
    return(NULL)
  }
  
  #grouping1_up$gene change : to -
  target_peaks <- gsub(":", "-", target_peaks)
  
  
  enriched.motifs <- FindMotifs(
    object = atac_object,
    features = target_peaks
  )
  enriched.motifs
  
  closest_features <- ClosestFeature(
    object = atac_object,
    regions = target_peaks
  )
  
  
  
  output_file = file.path(output_dir, sprintf("enrichment_results_%s_%s.rds", grouping_name, direction))
  # Save results
  saveRDS(enriched.motifs, output_file)
  write.csv(enriched.motifs, file.path(output_dir, sprintf("enrichment_results_%s_%s.csv", grouping_name, direction)), row.names = FALSE)
  
  #save closest_features
  write.csv(closest_features, file.path(output_dir, sprintf("closest_features_%s_%s.csv", grouping_name, direction)), row.names = FALSE)
  
  return(enriched.motifs)
}








# Create output directory
output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_analysis_results/sex/"
dir.create(output_dir, showWarnings = FALSE)

# Get unique groupings from sex_specific_regions$...
groupings <- unique(meta_data$assignments)
groupings <- c("Lactotrophs","Gonadotrophs","Stem_cells","Somatotrophs","Thyrotrophs","Melanotrophs","Corticotrophs")
directions <- c("male", "female")

# Initialize list to store all results
all_enrichment_results <- list()

# Loop through each grouping and direction
for(grouping in groupings) {
  for(direction in directions) {
    message(sprintf("Processing %s (%s)...", grouping, direction))
    
    result <- perform_motif_analysis(sex_specific_regions,
                                     grouping, 
                                     direction, 
                                     output_dir)
    
    if(!is.null(result)) {
      all_enrichment_results[[paste(grouping, direction, sep="_")]] <- result
    }
  }
}

# Save complete results object
saveRDS(all_enrichment_results, file.path(output_dir, "all_motif_results.rds"))








# Create summary table
summary_table <- data.frame()

for(name in names(all_enrichment_results)) {
  result <- all_enrichment_results[[name]]
  #keep pvalue < 0.05
  #result <- result[result$p.adjust < 0.05,]
  result$analysis <- name
  #if there are at least 5 enriched motifs add to summary table
  if(nrow(result) > 5) {
    #append to summary table
    summary_table <- rbind(summary_table, result)
    
  }
}

#extract rows where motif.name is Ar or ESR1
summary_table_hormones <- summary_table[summary_table$motif.name %in% c("Ar","ESR1"),]
summary_table_hormones$sex <- gsub("_.*", "", summary_table_hormones$analysis)
summary_table_hormones$cell_type <- gsub(".*_", "", summary_table_hormones$analysis)
summary_table_hormones <- summary_table_hormones %>% mutate(label_position = ifelse(motif.name == "Ar", 5, -5))

#make a plot

# Extract rows where motif.name is Ar or ESR1
summary_table_hormones <- summary_table[summary_table$motif.name %in% c("Ar","ESR1"),]
#rename ESR1 to Esr1
summary_table_hormones$motif.name <- gsub("ESR1", "Esr1", summary_table_hormones$motif.name)

# Extract cell type from analysis column
summary_table_hormones$cell_type <- gsub("_.*", "", summary_table_hormones$analysis)

# Make ESR1 fold enrichment negative for visualization (left side)
summary_table_hormones$fold.enrichment <- ifelse(
  summary_table_hormones$motif.name == "Esr1",
  -summary_table_hormones$fold.enrichment,
  summary_table_hormones$fold.enrichment
)

# Add label position for the motif names
summary_table_hormones <- summary_table_hormones %>% 
  mutate(label_position = ifelse(motif.name == "Ar", 
                                 fold.enrichment + max(fold.enrichment) * 0.1, 
                                 fold.enrichment - max(abs(fold.enrichment)) * 0.1))

# Calculate total enrichment (absolute values) for each cell type for sorting
cell_type_totals <- summary_table_hormones %>%
  group_by(cell_type) %>%
  summarize(total_enrichment = sum(abs(fold.enrichment))) %>%
  arrange(desc(total_enrichment))

# Reorder cell types by total enrichment
summary_table_hormones$cell_type <- factor(
  summary_table_hormones$cell_type,
  levels = cell_type_totals$cell_type
)

#summary_table_hormones
#only keep where pvalue <.05
summary_table_hormones <- summary_table_hormones[summary_table_hormones$pvalue < 0.05,]

# Make a plot
plt <- ggplot(summary_table_hormones, aes(x = cell_type, y = fold.enrichment, fill = motif.name)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = motif.name, y = label_position),
            color = "black",
            size = 5) +
  scale_fill_manual(values = c("Esr1" = "#FFA500", "Ar" = "#63B3ED")) +
  scale_y_continuous(
    labels = abs,
    limits = c(-max(abs(summary_table_hormones$fold.enrichment)) * 1.2, 
               max(abs(summary_table_hormones$fold.enrichment)) * 1.2)
  ) +
  coord_flip() +
  labs(
    title = "Estrogen and Androgen Receptor Motif Enrichment",
    x = "Cell Type",
    y = "Fold Enrichment",
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

# Save the plot
fold <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/"
ggsave(plt, filename = paste0(fold, "hormone_receptor_enrichment.png"), width = 4.5, height = 5.5, dpi = 300)
ggsave(plt, filename = paste0(fold, "hormone_receptor_enrichment.svg"), width = 4.5, height = 5.5, dpi = 300)

# Display the plot
plt


# Save summary table
write.csv(summary_table, 
          file.path(output_dir, "motif_analysis_summary.csv"), 
          row.names = FALSE)





library(tidyverse)
library(patchwork)
library(ggplot2)

# Set the output directory
output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_analysis_results/sex/"

# Load the Gonadotroph male and female motif lists
gonadotroph_male <- readRDS(file.path(output_dir, "enrichment_results_Gonadotrophs_male.rds"))
gonadotroph_female <- readRDS(file.path(output_dir, "enrichment_results_Gonadotrophs_female.rds"))

# Add a 'sex' column to each data frame
gonadotroph_male$sex <- "male"
gonadotroph_female$sex <- "female"

# Keep only necessary columns before merging
gonadotroph_male <- gonadotroph_male %>%
  dplyr::select(motif.name, pvalue, p.adjust, sex, fold.enrichment)
gonadotroph_female <- gonadotroph_female %>%
  dplyr::select(motif.name, pvalue, p.adjust, sex, fold.enrichment)

# Merge the two data frames
gonadotrophs_merged <- rbind(gonadotroph_male, gonadotroph_female)
#first name is chars before . or :
gonadotrophs_merged$first_name =  gsub("[:.-].*$", "", gonadotrophs_merged$motif.name)
#turn to CAPS
gonadotrophs_merged$first_name = toupper(gonadotrophs_merged$first_name)
#remove duplicates based on first_name keeping the one with the lowest p.adjust
gonadotrophs_merged <- gonadotrophs_merged %>%
  group_by(first_name) %>%
  #keep top 1
  slice_min(order_by = p.adjust, n = 1, with_ties = FALSE
  ) %>%
  ungroup()


# Calculate -log10 p.adjust and filter for p.adjust < 0.05
gonadotrophs_merged <- gonadotrophs_merged %>%
  mutate(p.adjust = p.adjust + 10**-300) %>%
  
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(neg_log10_p_adjust = -log10(p.adjust))

# Find the top 20 motifs based on the maximum -log10 p.adjust value between male and female
top_20_motifs <- gonadotrophs_merged %>%
  group_by(motif.name) %>%
  
  summarize(max_fold.enrichment = max(fold.enrichment, na.rm = TRUE)) %>%
  arrange(desc(max_fold.enrichment )) %>%
  slice_head(n = 30)

# Filter the merged data frame to keep only the top 20 motifs
plot_data <- gonadotrophs_merged %>%
  dplyr::filter(motif.name %in% c(top_20_motifs$motif.name))

# Reorder the motif.name factor based on the average -log10 p.adjust value
plot_data$motif.name <- factor(
  plot_data$motif.name,
  levels = unique(top_20_motifs$motif.name)
)

#if motif.name is NA fill with first_name
motif_names = as.character(plot_data$motif.name)
plot_data$motif.name <- ifelse(is.na(motif_names), plot_data$first_name,motif_names)

#sort levels by fold.enrichment
plot_data <- plot_data %>%
  arrange(desc(fold.enrichment)) %>%
  mutate(motif.name = factor(motif.name, levels = unique(motif.name)))

# Create the plot with male and female bars going in the same direction
plt_gonadotrophs <- ggplot(plot_data, aes(x = fold.enrichment, y = motif.name, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("male" = "#63B3ED", "female" = "#FFA500")) +
  labs(
    title = "Top 20 Enriched Motifs in Gonadotrophs",
    x = "Fold enrichment",
    y = "Motif Name",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

plt_gonadotrophs 

# Save the plot
ggsave(plt_gonadotrophs, filename = file.path(output_dir, "gonadotroph_top20_motifs.png"),
       width = 3.5, height = 5.5, dpi = 300)
ggsave(plt_gonadotrophs, filename = file.path(output_dir, "gonadotroph_top20_motifs.svg"),
       width = 3.5, height = 5.5, dpi = 300)

# Display the plot
print(plt_gonadotrophs)




library(tidyverse)
library(patchwork)
library(ggplot2)

# Set the output directory
output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_analysis_results/sex/"

# Load the Gonadotroph male and female motif lists
lactotroph_male<- readRDS(file.path(output_dir, "enrichment_results_Lactotrophs_male.rds"))
lactotroph_female <- readRDS(file.path(output_dir, "enrichment_results_Lactotrophs_female.rds"))

# Add a 'sex' column to each data frame
lactotroph_male$sex <- "male"
lactotroph_female$sex <- "female"

# Keep only necessary columns before merging
lactotroph_male <- lactotroph_male %>%
  dplyr::select(motif.name, pvalue, p.adjust, sex, fold.enrichment)
lactotroph_female <- lactotroph_female %>%
  dplyr::select(motif.name, pvalue, p.adjust, sex, fold.enrichment)

# Merge the two data frames
lactotrophs_merged <- rbind(lactotroph_male, lactotroph_female)
#first name is chars before . or :
lactotrophs_merged$first_name =  gsub("[:.-].*$", "", lactotrophs_merged$motif.name)
#turn to CAPS
lactotrophs_merged$first_name = toupper(lactotrophs_merged$first_name)
#remove duplicates based on first_name keeping the one with the lowest p.adjust
lactotrophs_merged <- lactotrophs_merged %>%
  group_by(first_name) %>%
  #keep top 1
  slice_min(order_by = p.adjust, n = 1, with_ties = FALSE
            ) %>%
  ungroup()
  

# Calculate -log10 p.adjust and filter for p.adjust < 0.05
lactotrophs_merged <- lactotrophs_merged %>%
  mutate(p.adjust = p.adjust + 10**-300) %>%
  
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(neg_log10_p_adjust = -log10(p.adjust))

# Find the top 20 motifs based on the maximum -log10 p.adjust value between male and female
top_20_motifs <- lactotrophs_merged %>%
  group_by(motif.name) %>%
  
  summarize(max_fold.enrichment = max(fold.enrichment, na.rm = TRUE)) %>%
  arrange(desc(max_fold.enrichment )) %>%
  slice_head(n = 30)

# Filter the merged data frame to keep only the top 20 motifs
plot_data <- lactotrophs_merged %>%
  dplyr::filter(motif.name %in% c(top_20_motifs$motif.name))

# Reorder the motif.name factor based on the average -log10 p.adjust value
plot_data$motif.name <- factor(
  plot_data$motif.name,
  levels = unique(top_20_motifs$motif.name)
)

#if motif.name is NA fill with first_name
motif_names = as.character(plot_data$motif.name)
plot_data$motif.name <- ifelse(is.na(motif_names), plot_data$first_name,motif_names)

#sort levels by fold.enrichment
plot_data <- plot_data %>%
  arrange(desc(fold.enrichment)) %>%
  mutate(motif.name = factor(motif.name, levels = unique(motif.name)))

# Create the plot with male and female bars going in the same direction
plt_lactotrophs <- ggplot(plot_data, aes(x = fold.enrichment, y = motif.name, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("male" = "#63B3ED", "female" = "#FFA500")) +
  labs(
    title = "Top 20 Enriched Motifs in Lactotrophs",
    x = "Fold enrichment",
    y = "Motif Name",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

plt_lactotrophs 

# Save the plot
ggsave(plt_lactotrophs, filename = file.path(output_dir, "lactotroph_top20_motifs.png"),
       width = 3.5, height = 5.5, dpi = 300)
ggsave(plt_lactotrophs, filename = file.path(output_dir, "lactotroph_top20_motifs.svg"),
       width = 3.5, height = 5.5, dpi = 300)

# Display the plot
print(plt_lactotrophs)




library(tidyverse)
library(patchwork)
library(ggplot2)

# Set the output directory
output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/motif_analysis_results/sex/"

# Load the Gonadotroph male and female motif lists
sc_male <- readRDS(file.path(output_dir, "enrichment_results_Stem_cells_male.rds"))
sc_female <- readRDS(file.path(output_dir, "enrichment_results_Stem_cells_female.rds"))

# Add a 'sex' column to each data frame
sc_male$sex <- "male"
sc_female$sex <- "female"

# Keep only necessary columns before merging
sc_male <- sc_male %>%
  dplyr::select(motif.name, pvalue, p.adjust, sex, fold.enrichment)
sc_female <- sc_female %>%
  dplyr::select(motif.name, pvalue, p.adjust, sex, fold.enrichment)

# Merge the two data frames
scs_merged <- rbind(sc_male, sc_female)
#first name is chars before . or :
scs_merged$first_name =  gsub("[:.-].*$", "", scs_merged$motif.name)
#turn to CAPS
scs_merged$first_name = toupper(scs_merged$first_name)
#remove duplicates based on first_name keeping the one with the lowest p.adjust
scs_merged <- scs_merged %>%
  group_by(first_name) %>%
  #keep top 1
  slice_min(order_by = p.adjust, n = 1, with_ties = FALSE
  ) %>%
  ungroup()


# Calculate -log10 p.adjust and filter for p.adjust < 0.05
scs_merged <- scs_merged %>%
  mutate(p.adjust = p.adjust + 10**-300) %>%
  
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(neg_log10_p_adjust = -log10(p.adjust))

# Find the top 20 motifs based on the maximum -log10 p.adjust value between male and female
top_20_motifs <- scs_merged %>%
  group_by(motif.name) %>%
  
  summarize(max_fold.enrichment = max(fold.enrichment, na.rm = TRUE)) %>%
  arrange(desc(max_fold.enrichment )) %>%
  slice_head(n = 30)

# Filter the merged data frame to keep only the top 20 motifs
plot_data <- scs_merged %>%
  dplyr::filter(motif.name %in% c(top_20_motifs$motif.name))

# Reorder the motif.name factor based on the average -log10 p.adjust value
plot_data$motif.name <- factor(
  plot_data$motif.name,
  levels = unique(top_20_motifs$motif.name)
)

#if motif.name is NA fill with first_name
motif_names = as.character(plot_data$motif.name)
plot_data$motif.name <- ifelse(is.na(motif_names), plot_data$first_name,motif_names)

#sort levels by fold.enrichment
plot_data <- plot_data %>%
  arrange(desc(fold.enrichment)) %>%
  mutate(motif.name = factor(motif.name, levels = unique(motif.name)))

# Create the plot with male and female bars going in the same direction
plt_scs <- ggplot(plot_data, aes(x = fold.enrichment, y = motif.name, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("male" = "#63B3ED", "female" = "#FFA500")) +
  labs(
    title = "Top 20 Enriched Motifs in scs",
    x = "Fold enrichment",
    y = "Motif Name",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

plt_scs 

# Save the plot
ggsave(plt_scs, filename = file.path(output_dir, "sc_top20_motifs.png"),
       width = 3.5, height = 5.5, dpi = 300)
ggsave(plt_scs, filename = file.path(output_dir, "sc_top20_motifs.svg"),
       width = 3.5, height = 5.5, dpi = 300)

# Display the plot
print(plt_scs)




