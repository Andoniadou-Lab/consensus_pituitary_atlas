
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


# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation


L <- makeContrastsDream(form, meta_data,
                        contrasts = c( 
                          
                          #Stem cells
                          Stem_cells_age = "assignmentsStem_cells.Age_numeric",
                          Corticotrophs = "assignmentsCorticotrophs.Age_numeric",
                          Melanotrophs = "assignmentsMelanotrophs.Age_numeric",
                          Thyrotrophs = "assignmentsThyrotrophs.Age_numeric",
                          Somatotrophs = "assignmentsSomatotrophs.Age_numeric",
                          Lactotrophs = "assignmentsLactotrophs.Age_numeric",
                          assignmentsPituicytes = "assignmentsPituicytes.Age_numeric",
                          assignmentsImmune_cells = "assignmentsImmune_cells.Age_numeric",
                          assignmentsMesenchymal_cells = "assignmentsMesenchymal_cells.Age_numeric",
                          assignmentsEndothelial_cells = "assignmentsEndothelial_cells.Age_numeric"
                          
                        )
)


# fit dream model with contrasts
fit <- dream(vobjDream, form, meta_data, L)

head(coef(fit))

fit <- dream(vobjDream, form, meta_data)
fit <- eBayes(fit)




background_genes <- rownames(dge)
#save to /Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/enrichment/aging_enrichment/aging_analysis_summary.csv
write.csv(background_genes, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/enrichment/aging_enrichment/sc_aging_background_genes.csv", row.names = FALSE)

#####
#Age effect
#####

contrast <- makeContrasts(
  #assignmentsStem_cells - assignmentsLactotrophs, 
  #assignmentsGonadotrophs.Age_numeric,
  assignmentsStem_cells.Age_numeric,
  levels = design
)

# Compute differential expression
diff_expression <- contrasts.fit(fit, contrast)
diff_expression <- eBayes(diff_expression, robust=TRUE)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf, adjust.method="bonf", lfc=0.5)
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))

#how many sig
print(nrow(top_genes %>% filter(adj.P.Val < 0.05 )))


cpdb <- read_csv("/Users/k23030440/epitome_code/epitome/data/gene_group_annotation/v_0.02//cpdb.csv")
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

cpdb <- read_csv("/Users/k23030440/epitome_code/epitome/data/gene_group_annotation/v_0.02//cpdb.csv")
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
  top_genes <- topTable(diff_expression, number = Inf,adjust.method="bonf")
  top_genes$genes <- rownames(top_genes)
  #print cell type
  print(get_cell_type(coef))
  #print number of significant genes
  print(nrow(top_genes %>% filter(adj.P.Val < 0.05 )))
  
  # Add cell type column
  top_genes$cell_type <- get_cell_type(coef)
  
  # Store results
  all_results[[coef]] <- top_genes
}

# Combine all results into single dataframe
final_results <- do.call(rbind, all_results)

# Reset row names
rownames(final_results) <- NULL



#save to /Users/k23030440/epitome_code/epitome/aging/v_0.02
version = "v_0.02"
write.csv(final_results, paste0("/Users/k23030440/epitome_code/epitome/data/aging/v_0.02/aging_genes.csv"))


final_results <- read.csv(paste0("/Users/k23030440/epitome_code/epitome/data/aging/v_0.02/aging_genes.csv"), header = TRUE, stringsAsFactors = FALSE)

#keep only padj < 0.05
final_results_pan_pit <- final_results[final_results$adj.P.Val < 0.05,]
#filter log2fc for abs 0.5
final_results_pan_pit <- final_results_pan_pit[abs(final_results_pan_pit$logFC) >= 0.5,]
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

#save these genes as csv
#/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/DE/pan_pituitary_ageing.csv
write.csv(final_results_pan_pit, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/pan_pituitary_ageing.csv")
#clipboard gene names
unique_pan_pit_genes = final_results_pan_pit$genes
library(clipr)
write_clip(paste(unique_pan_pit_genes, collapse = "\n"))

#I have this df called final_results with columns logFC, genes, adj.P.Val and cell_type. Make a heatmap, where tiles are outline wiht red if significant and colored with light to dark blue based on logFC. Make columns cell_types and rows genes. Specifically look at Ifit3, Mki67 and Lef1


# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)

# Highlight specific genes of interest
genes_of_interest <- c("Wnt10a","Wnt16","Wnt6","Fzd2","Fzd9","Fzd10","Wif1","Lef1","Rspo4","Sfrp1","Sfrp2","Sfrp5","Apcdd1","Tnfrsf19","Sostdc1",
                       "Cdk1","Cdk2","Ccna1","Ccnb1","H19","E2f8","E2f1","Mki67","Top2a","Ube2c","Aurkb",
                       "Lcn2","Cxcl13","Ccl5","H2-Aa","H2-K1","H2-Ab1","Ifit3","Ifit3b","C1qc","C1qa","Il1b","Il6","Il18","Il17rc","Il15ra","Igha")

#keep only where genes in final_results_pan_pit
genes_of_interest <- genes_of_interest[genes_of_interest %in% final_results_pan_pit$genes]
# Create the heatmap
plot <- final_results %>%
  # Filter to only include genes of interest
  filter(genes %in% genes_of_interest) %>%
  # Create a significance flag and preserve input order
  dplyr::mutate(is_significant = adj.P.Val < 0.05,
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
                       linewidth = 2
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
#plot_data <- plot_data[!grepl("Immune_cells", plot_data$cell_type), ]
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
  dplyr::mutate(total_diff_exp = sum(abs(count))) %>% # Total absolute differentially expressed genes
  ungroup() %>%
  arrange(desc(total_diff_exp)) %>%            # Order by total genes
  dplyr::mutate(cell_type = factor(cell_type, levels = unique(cell_type))) # Maintain order

# Add a column for fixed label positions
plot_data_long <- plot_data_long %>%
  dplyr::mutate(label_position = ifelse(count > 0, 400, -400)) # Labels at fixed positions

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
  width = 5.5,
  height = 6,
  dpi = 300
)

# SVG format
ggsave(
  filename = file.path(output_dir, "age_related_genes.svg"),
  plot = plt,
  width = 5.5,
  height = 6
)

print(plt)




#############
#############
# Age related patterns
#############
#############

#load ageing results /Users/k23030440/epitome_code/epitome/data/aging/v_0.01/aging_genes.csv
aging = "/Users/k23030440/epitome_code/epitome/data/aging/v_0.01/aging_genes.csv"
aging = read.csv(aging)
aging

#do model selection on steady, change-steady, steady-change
#load normalised counts
Matrix::writeMM(normalized_expression, paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/normalized_data.mtx"))
write_csv(meta_data_original, paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/meta_data.csv"))
write.table(rownames(normalized_expression), paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(normalized_expression), paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/datasets.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

rownames(normalized_expression)

#LOADING
library(Matrix)
library(readr)

version="v_0.01"
# 1. Load the matrix
normalized_expression <- Matrix::readMM(paste0("/Users/k23030440/epitome_code/epitome/data/expression/", version, "/normalized_data.mtx"))

# 2. Load the metadata
meta_data <- read_csv(paste0("/Users/k23030440/epitome_code/epitome/data/expression/", version, "/meta_data.csv"))

# 3. Load and assign Genes (Row names)
genes <- read.table(paste0("/Users/k23030440/epitome_code/epitome/data/expression/", version, "/genes.txt"), stringsAsFactors = FALSE)$V1
rownames(normalized_expression) <- genes

# 4. Load and assign Cells/Datasets (Column names)
datasets <- read.table(paste0("/Users/k23030440/epitome_code/epitome/data/expression/", version, "/datasets.txt"), stringsAsFactors = FALSE)$V1
colnames(normalized_expression) <- datasets

version="v_0.02"


#keep only where Modality is "sc"
normalized_expression <- normalized_expression[, meta_data$Modality == "sc"]
meta_data <- meta_data[meta_data$Modality == "sc", ]
# Check the result
print(dim(normalized_expression)) 



#define a function get values, that takes cell_type and gene_name and returns two vectors: age and gene exp
get_values <- function(cell_type, gene_name, meta, expr_matrix) {
  # Get indices for the specific cell type
  idx <- which(meta$assignments == cell_type)
  
  # Return age and expression as a list
  return(list(
    age = log10(meta$Age_numeric[idx]), # Log-transform age as requested
    expression = log10(as.numeric(expr_matrix[gene_name, idx])+1)
  ))
}



data_list <- get_values("Stem_cells", "H19", meta_data, normalized_expression)

analyze_gene_pattern <- function(data_list, verbose = TRUE) {
  library(segmented)
  
  # 1. Prepare Data
  df <- data.frame(age = data_list$age, expression = data_list$expression)
  df <- df[order(df$age), ]
  
  # 2. Fit Models
  fit_linear <- lm(expression ~ age, data = df)
  fit_seg <- try(segmented(fit_linear, seg.Z = ~age, 
                           psi = median(df$age)), silent = TRUE)
  
  # 3. Model Selection Logic
  # We check if segmented fit exists and if it improves AIC by at least 2
  if (!inherits(fit_seg, "try-error") && (AIC(fit_linear) - AIC(fit_seg) > 2)) {
    
    # Extract slopes
    slopes <- slope(fit_seg)$age
    s1 <- slopes[1, "Est."]
    s2 <- slopes[2, "Est."]
    
    # Extract breakpoint (psi)
    # fit_seg$psi is a matrix; [1, 2] is the estimated breakpoint value
    breakpoint <- fit_seg$psi[1, 2]
    
    # Categorize
    if (abs(s1) > abs(s2)) {
      final_scenario <- ifelse(s1 > 0, "increase-stagnate", "decrease-stagnate")
    } else {
      final_scenario <- ifelse(s2 > 0, "stagnate-increase", "stagnate-decrease")
    }
    
  } else {
    # Linear model wins
    lin_slope <- coef(fit_linear)[2]
    final_scenario <- ifelse(lin_slope > 0, "increase steadily", "decrease steadily")
    breakpoint <- -10
  }
  
  # 4. Output
  
  if  (verbose) {
    print(paste("The scenario for this gene is:", final_scenario))
    if(breakpoint != -10) {
      print(paste("Breakpoint detected at log age:", round(breakpoint, 3)))
    }
  }
  
  return(list(
    scenario = final_scenario,
    breakpoint = breakpoint
  ))
}


plot_gene_models <- function(data_list, gene_label = "Gene") {
  library(segmented)
  
  # 1. Prepare and Sort Data
  df <- data.frame(age = data_list$age, expression = data_list$expression)
  df <- df[order(df$age), ]
  
  # 2. Fit the Models
  fit_steady <- lm(expression ~ age, data = df)
  fit_seg <- try(segmented(fit_steady, seg.Z = ~age, 
                           psi = median(df$age)), silent = TRUE)
  
  # 3. Create the Base Plot
  plot(df$age, df$expression, 
       pch = 16, 
       col = rgb(0.2, 0.2, 0.2, 0.3), 
       main = paste("Competing Models for", gene_label),
       xlab = "Log Age", 
       ylab = "Normalized Expression")
  
  # 4. Add Steady Model (Blue Dashed)
  abline(fit_steady, col = "royalblue", lwd = 2, lty = 2)
  
  # 5. Add Segmented Model (Red Solid)
  if (!inherits(fit_seg, "try-error")) {
    plot(fit_seg, add = TRUE, col = "firebrick", lwd = 3)
    
    # Add a vertical line at the detected breakpoint
    # fit_seg$psi[1, 2] is the coordinate on the x-axis
    abline(v = fit_seg$psi[1, 2], col = "darkgrey", lty = 3)
    
    # Determine type for plot sub-label
    slopes <- slope(fit_seg)$age
    type_label <- if(abs(slopes[1,1]) > abs(slopes[2,1])) "Early Change" else "Late Change"
    mtext(paste("Pattern Type:", type_label), side = 3, line = 0.5, cex = 0.8, col = "firebrick")
  }
  
  # 6. Add Legend
  legend("topleft", 
         legend = c("Steady Linear", "Segmented", "Breakpoint"),
         col = c("royalblue", "firebrick", "darkgrey"), 
         lty = c(2, 1, 3), 
         lwd = 2,
         bg = "white")
}


# 2. Analyze
gene="Gh"
cell_type = "Stem_cells"

dim(aging)
#keep onyl where adj.P.Val < 0.01
aging = aging[aging$adj.P.Val < 0.05, ]
#keep onyl Stem_cells for now
aging = aging[aging$cell_type == cell_type, ]
dim(aging)

data_h19 <- get_values(cell_type,gene, meta_data, normalized_expression)
result <- analyze_gene_pattern(data_h19)
print(result$scenario)
print(result$breakpoint)
plot_gene_models(data_h19, gene_label = gene)



# 1. Initialize an empty list to store results
results_list <- list()

# 2. Get unique cell types from the aging dataframe
unique_cell_types <- unique(aging$cell_type)

for (ct in unique_cell_types) {
  message(paste("Processing cell type:", ct))
  
  # Filter aging DF for significant genes in this cell type
  sig_genes_df <- aging[aging$cell_type == ct,]
  genes_to_test <- sig_genes_df$genes
  
  for (g in genes_to_test) {
    # Wrapped in try to prevent the whole loop from crashing if one gene fails
    res <- try({
      # Get data
      d_list <- get_values(ct, g, meta_data, normalized_expression)
      
      # Analyze pattern (using the function we built)
      analysis <- analyze_gene_pattern(d_list, verbose = FALSE)
      
      # Store in a temporary data frame
      data.frame(
        cell_type = ct,
        gene = g,
        change_pattern = analysis$scenario,
        breakpoint = analysis$breakpoint,
        stringsAsFactors = FALSE
      )
    }, silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      results_list[[length(results_list) + 1]] <- res
    }
  }
}

# 3. Combine all results into one Data Frame
final_summary_df <- do.call(rbind, results_list)

# 4. View results
head(final_summary_df)








#########
###
##

###
meta_data <- as.data.frame(adata$obs)
#keep only unique SRA_IDs
meta_data <- meta_data %>% distinct(SRA_ID, .keep_all = TRUE)
# Mann-Whitney U Test
wilcox_result <- wilcox.test(Age_numeric ~ Comp_sex, data = meta_data)

# Print results
print(wilcox_result)

# Quick visualization to check the overlap
boxplot(Age_numeric ~ Comp_sex, data = meta_data, 
        main = "Age Distribution by Sex", 
        ylab = "Age", col = c("lightblue", "lightpink"))

#repeat on log10 scale
meta_data$log10_Age <- log10(meta_data$Age_numeric + 1) # Adding 1 to avoid log(0)
wilcox_result_log <- wilcox.test(log10_Age ~ Comp_sex, data = meta_data)
boxplot(log10_Age ~ Comp_sex, data = meta_data, 
        main = "Log10 Age Distribution by Sex", 
        ylab = "Log10 Age", col = c("lightblue", "lightpink"))
print(wilcox_result_log)





# 1. Ensure Age is integer
meta_data$Age_numeric <- as.integer(meta_data$Age_numeric)

meta_data$Age_bins <- cut(meta_data$Age_numeric, 
                          breaks = c(0, 100, 200,2000), 
                          labels = c("0-100", "101-200", "201+"),
                          right = FALSE)

# 3. Create a contingency table (Removing empty levels)
age_sex_table <- table(meta_data$Age_bins, meta_data$Comp_sex)
# Drop bins that have zero observations across both sexes
age_sex_table <- age_sex_table[rowSums(age_sex_table) > 0, ]

# 4. Perform Chi-Square Test
chisq_result <- chisq.test(age_sex_table)

# 5. If you still get a warning, use Fisher's Exact Test
# This is more accurate for small counts or unbalanced tables
fisher_result <- fisher.test(age_sex_table, simulate.p.value = TRUE)

# Print results
print(age_sex_table)
print(chisq_result)
print(fisher_result)


