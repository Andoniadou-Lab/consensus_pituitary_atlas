#Loading dependencies
library(limma)
library(edgeR)
library(readr)
library(tibble) 
library(dplyr)
library(Matrix)

# --- Set Working Directory ---
setwd("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/bulk validations/folder_with_bulk_results/")

figs_folder <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/bulk validations/"

# --- Load Metadata ---
metadata <- read.csv("metadata.csv")
# --- Initialize an Empty List to Store Matrices ---
counts_list <- list()

# --- Loop and Load Data ---
for (sra_id in metadata$SRA_ID) {
  data_dir <- file.path(sra_id)
  matrix_file <- file.path(data_dir, "cells_x_genes.total.mtx")
  genes_file <- file.path(data_dir, "cells_x_genes.genes.names.txt")
  
  # --- Load Data with File Checks ---
  if (file.exists(matrix_file) && file.exists(genes_file)) {
    counts <- readMM(matrix_file)
    #transpose
    counts <- t(counts)
    gene_names <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)$V1
    
    # Dimension Check
    if (nrow(counts) != length(gene_names)) {
      stop(paste("Number of rows in matrix and length of gene names do not match for", sra_id))
    }
    
    # Make sure the counts are a single column.
    if (ncol(counts) != 1) {
      stop(paste("Matrix for", sra_id, "should have only 1 column, but has", ncol(counts)))
    }
    rownames(counts) <- gene_names
    colnames(counts) <- sra_id  # Use SRA_ID as the column name (sample name)
    
    counts_list[[sra_id]] <- counts  # Add the matrix to the list
    
    
  } else {
    warning(paste("Matrix or genes file not found for", sra_id, ". Skipping."))
  }
}

# --- Combining Matrices ---
# Get all unique genes.
all_genes <- unique(unlist(lapply(counts_list, rownames)))

# Create an empty matrix with all genes and sample IDs as columns
combined_counts <- matrix(0, nrow = length(all_genes), ncol = length(counts_list))
rownames(combined_counts) <- all_genes
colnames(combined_counts) <- names(counts_list)

# Fill in the combined matrix.
for (sra_id in names(counts_list)) {
  counts_matrix <- counts_list[[sra_id]]
  combined_counts[rownames(counts_matrix), sra_id] <- counts_matrix[, 1] # Use [,1] to handle single-column matrices.
}

# --- Create DGEList Object ---
dge <- DGEList(counts = combined_counts, samples = metadata)
dim(dge)
#keep where metadata$Author is Hou et al., (2022)
dge_hou <- dge[, metadata$Author == "Hou et al., (2022)"] #whole pituitaries
#keep where Duncan et al. (2023)
dge_duncan <- dge[, metadata$Author == "Duncan et al. (2023)"] #Corticotrophs
#keep where Alonso et al., (2023)
dge_alonso <- dge[, metadata$Author == "Alonso et al., (2023)"] #Hpg mice

#in dge_alonso keep only where Conditions is placebo
dge_alonso <- dge[, metadata$Author == "Alonso et al., (2023)" & metadata$Conditions == "placebo"]

#and conditions doesnt start with Nr5a1-PE
dge_shima <- dge[, metadata$Author == "Shima et al. (2022)" & grepl("Nr5a1-PE", metadata$Conditions) == FALSE] #Gonadotrophs


#in metadata remove where its alonso and not placebo
metadata <- metadata[!(metadata$Author == "Alonso et al., (2023)" & metadata$Conditions != "placebo"), ]
metadata <- metadata[!(metadata$Author == "Shima et al. (2022)" & grepl("Nr5a1-PE", metadata$Conditions)), ] 

#calcnormfactors for all samples
dge_hou <- calcNormFactors(dge_hou)
dge_duncan <- calcNormFactors(dge_duncan)
dge_alonso <- calcNormFactors(dge_alonso)
dge_shima <- calcNormFactors(dge_shima)


# Create the design matrix for all samples
design_hou <- model.matrix(~ 0 + Sex, data = metadata[metadata$Author == "Hou et al., (2022)", ])
design_duncan <- model.matrix(~ 0 + Sex, data = metadata[metadata$Author == "Duncan et al. (2023)", ])
design_alonso <- model.matrix(~ 0 + Sex, data = metadata[metadata$Author == "Alonso et al., (2023)", ])
design_shima <- model.matrix(~ 0 + Sex, data = metadata[metadata$Author == "Shima et al. (2022)", ])


#make syntactically correct names for all samples
colnames(design_hou) <- make.names(colnames(design_hou))
colnames(design_duncan) <- make.names(colnames(design_duncan))
colnames(design_alonso) <- make.names(colnames(design_alonso))
colnames(design_shima) <- make.names(colnames(design_shima))

#filterByExpr with default settings for all samples
keep_hou <- filterByExpr(dge_hou, design_hou)
keep_duncan <- filterByExpr(dge_duncan, design_duncan)
keep_alonso <- filterByExpr(dge_alonso, design_alonso)
keep_shima <- filterByExpr(dge_shima, design_shima)

#keeping only the chosen genes
dge_hou <- dge_hou[keep_hou, ]
dge_duncan <- dge_duncan[keep_duncan, ]
dge_alonso <- dge_alonso[keep_alonso, ]
dge_shima <- dge_shima[keep_shima, ]

#printing dimensions
print(dim(dge_hou))
print(dim(dge_duncan))
print(dim(dge_alonso))
print(dim(dge_shima))


#running voom on all datasets
fit_hou <- voom(dge_hou, design_hou, plot=TRUE)
fit_duncan <- voom(dge_duncan, design_duncan, plot=TRUE)
fit_alonso <- voom(dge_alonso, design_alonso, plot=TRUE)
fit_shima <- voom(dge_shima, design_shima, plot=TRUE)

#running lmFit on all datasets
fit_hou <- lmFit(fit_hou, design_hou)
fit_duncan <- lmFit(fit_duncan, design_duncan)
fit_alonso <- lmFit(fit_alonso, design_alonso)
fit_shima <- lmFit(fit_shima, design_shima)

#Make dataset specific contrasts
# Create contrast

contrast_hou <- makeContrasts(
  Sexfemale - Sexmale,
  levels = design_hou
)

contrast_duncan <- makeContrasts(
  Sexfemale - Sexmale,
  levels = design_duncan
)

contrast_alonso <- makeContrasts(
  Sexfemale - Sexmale,
  levels = design_alonso
)

contrast_shima <- makeContrasts(
  Sexfemale - Sexmale,
  levels = design_shima
)

# Compute differential expression
diff_hou <- contrasts.fit(fit_hou, contrast_hou)
diff_duncan <- contrasts.fit(fit_duncan, contrast_duncan)
diff_alonso <- contrasts.fit(fit_alonso, contrast_alonso)
diff_shima <- contrasts.fit(fit_shima, contrast_shima)

#running eBayes
diff_hou  <- eBayes(diff_hou )
diff_duncan  <- eBayes(diff_duncan )
diff_alonso  <- eBayes(diff_alonso )
diff_shima  <- eBayes(diff_shima )

# Get top differentially expressed genes
top_genes_hou <- topTable(diff_hou, coef = 1, number = Inf)
top_genes_duncan <- topTable(diff_duncan, coef = 1, number = Inf)
top_genes_alonso <- topTable(diff_alonso, coef = 1, number = Inf)
top_genes_shima <- topTable(diff_shima, coef = 1, number = Inf)

#*-1 the logFC for hou, duncan, alonso
top_genes_hou$logFC <- -1 * top_genes_hou$logFC
top_genes_duncan$logFC <- -1 * top_genes_duncan$logFC
top_genes_alonso$logFC <- -1 * top_genes_alonso$logFC
top_genes_shima$logFC <- -1 * top_genes_shima$logFC


# --- Save Results ---
saveRDS(fit_hou, "fit_hou.rds")
saveRDS(fit_duncan, "fit_duncan.rds")
saveRDS(fit_alonso, "fit_alonso.rds")
saveRDS(fit_shima, "fit_shima.rds")



library(ggrepel)
library(ggbreak)

# Defining the key plotting function for volcano plots that will display genes significant in the CPA, but plotted with bulk data
plot_consistent_volcano_with_gap <- function(df, y_gap_min = NULL, y_gap_max = NULL, show_ns = FALSE, p_val_threshold = 0.05,
                                             color1="blue", color2="purple", highlight_genes = NULL
) {
  consistent_sign_down <- rownames(df[(df$logFC < 0) &
                                        (df$gene_direction == "Down"),])
  consistent_sign_up <- rownames(df[(df$logFC > 0) &
                                      (df$gene_direction == "Up"),])
  
  # Calculate the number of consistent genes and the total in each direction
  n_consistent_down <- length(consistent_sign_down)
  total_down <- sum(df$gene_direction == "Down")
  
  n_consistent_up <- length(consistent_sign_up)
  total_up <- sum(df$gene_direction == "Up")
  
  print(paste0("Number of consistent down-regulated genes: ", n_consistent_down))
  print(paste0("Total down-regulated genes: ", total_down))
  print(paste0("Number of consistent up-regulated genes: ", n_consistent_up))
  print(paste0("Total up-regulated genes: ", total_up))
  
  # Calculate number of hits above p_val_threshold for each direction
  n_consistent_down_above_thresh <- sum((df$logFC < 0) & (df$gene_direction == "Down") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  total_down_above_thresh <- sum((df$gene_direction == "Down") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  
  n_consistent_up_above_thresh <- sum((df$logFC > 0) & (df$gene_direction == "Up") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  total_up_above_thresh <- sum((df$gene_direction == "Up") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  
  
  print(paste0("Number of consistent down-regulated genes above p-value threshold: ", n_consistent_down_above_thresh))
  print(paste0("Total down-regulated genes above p-value threshold: ", total_down_above_thresh))
  #for up
  print(paste0("Number of consistent up-regulated genes above p-value threshold: ", n_consistent_up_above_thresh))
  print(paste0("Total up-regulated genes above p-value threshold: ", total_up_above_thresh))
  
  # Calculate max y value before trimming for annotation placement
  max_y <- max(-log10(df$adj.P.Val), na.rm = TRUE)
  
  # Create the initial volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    # Plot "NS" points first (in the background) if show_ns is TRUE
    {if (show_ns) geom_point(data = subset(df, gene_direction == "NS"),
                             aes(color = gene_direction), size = 1.5, alpha = 0.5)} +
    # Plot "Up" and "Down" points on top
    geom_point(data = subset(df, gene_direction %in% c("Up", "Down")),
               aes(color = gene_direction), size = 1.5, alpha = 0.5) +
    scale_color_manual(values = c("Up" = color2, "Down" = color1, "NS" = "grey")) +
    labs(
      title = "Volcano Plot of DE Analysis",
      x = "log2 Fold Change (logFC)",
      y = "-log10(Adjusted P-value)"
    ) +
    theme_bw() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "black")+
    
    # Annotate for Down-regulated genes
    annotate("text", x = max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.3, hjust = 0.5,
             label = paste0(n_consistent_down, " (", n_consistent_down_above_thresh, ") consistent"), color = color2) +
    annotate("text", x = max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.2, hjust = 0.5,
             label = paste0("out of ", total_down, " (", total_down_above_thresh, ")"), color = color2) +
    
    # Annotate for Up-regulated genes
    annotate("text", x = -max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.3, hjust = 0.5,
             label = paste0(n_consistent_up, " (", n_consistent_up_above_thresh, ") consistent"), color = color1) +
    annotate("text", x = -max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.2, hjust = 0.5,
             label = paste0("out of ", total_up, " (", total_up_above_thresh, ")"), color =color1)
  
  # Highlight specified genes with a strong red dot
  if (!is.null(highlight_genes)) {
    highlight_data <- df[rownames(df) %in% highlight_genes, ]
    print(highlight_data)
    if (nrow(highlight_data) > 0) {
      p <- p + geom_point(data = highlight_data, aes(x = logFC, y = -log10(adj.P.Val)),
                          color = "red", size = 3)
    }
  }
  
  # Add y-axis break if y_gap_min and y_gap_max are provided
  if (!is.null(y_gap_min) && !is.null(y_gap_max)) {
    p <- p + scale_y_break(c(y_gap_min, y_gap_max),scales = 0.5)
  }
  
  return(p)
}


##############
####Duncan Corticotroph validations!
##############

#load /Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv
sexually_dimorphic_genes <- read.csv("/Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv")
sexually_dimorphic_genes_corticotrophs <- sexually_dimorphic_genes[sexually_dimorphic_genes$cell_type == "Corticotrophs",]
up_genes <- sexually_dimorphic_genes_corticotrophs[sexually_dimorphic_genes_corticotrophs$logFC > 0,]$gene
down_genes <- sexually_dimorphic_genes_corticotrophs[sexually_dimorphic_genes_corticotrophs$logFC < 0,]$gene


library(ggplot2)
library(dplyr)

# Convert top_genes_duncan to a data frame if it's not already
top_genes_duncan_df <- as.data.frame(top_genes_duncan)

# Add a column to indicate gene direction (up, down, or NS - not significant/not in our lists)
top_genes_duncan_df <- top_genes_duncan_df %>%
  mutate(gene_direction = case_when(
    rownames(.) %in% up_genes ~ "Up",
    rownames(.) %in% down_genes  ~ "Down",
    TRUE ~ "NS"
  ))

volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_duncan_df,show_ns= FALSE,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)



#save plot to figs_folder
ggsave(paste0(figs_folder, "duncan_volcano_plot_sex_de_min_3.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "duncan_volcano_plot_sex_de_min_3.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)



volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_duncan_df,show_ns= TRUE,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "duncan_volcano_plot_sex_de_min_3_NS.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "duncan_volcano_plot_sex_de_min_3_NS.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)




##############
####Shima Gonadotroph validations!
##############

#load /Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv
sexually_dimorphic_genes <- read.csv("/Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv")
sexually_dimorphic_genes_corticotrophs <- sexually_dimorphic_genes[sexually_dimorphic_genes$cell_type == "Gonadotrophs",]
up_genes <- sexually_dimorphic_genes_corticotrophs[sexually_dimorphic_genes_corticotrophs$logFC > 0,]$gene
down_genes <- sexually_dimorphic_genes_corticotrophs[sexually_dimorphic_genes_corticotrophs$logFC < 0,]$gene


library(ggplot2)
library(dplyr)

# Convert top_genes_duncan to a data frame if it's not already
top_genes_shima_df <- as.data.frame(top_genes_shima)

# Add a column to indicate gene direction (up, down, or NS - not significant/not in our lists)
top_genes_shima_df <- top_genes_shima_df %>%
  mutate(gene_direction = case_when(
    rownames(.) %in% up_genes ~ "Up",
    rownames(.) %in% down_genes  ~ "Down",
    TRUE ~ "NS"
  ))



volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_shima_df,show_ns= FALSE,color1="#ffa500",color2="#63b3ed",highlight_gene = c("Ppp1rc1","Gpr101","Grem1","Fshb"))
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "shima_volcano_plot_sex_de_min_3.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "shima_volcano_plot_sex_de_min_3.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)



volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_shima_df,show_ns= TRUE,color1="#ffa500",color2="#63b3ed",highlight_gene = c("Fshb"))
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "shima_volcano_plot_sex_de_min_3_NS.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "shima_volcano_plot_sex_de_min_3_NS.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)


##############
####Alonso (hpg mice) pituitary-wide validations!
##############


#load /Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv
sexually_dimorphic_genes <- read.csv("/Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv")
sexually_dimorphic_genes_pit_wide <- sexually_dimorphic_genes[sexually_dimorphic_genes$occurs > 1,]
sexually_dimorphic_genes_pit_wide


up_genes <- sexually_dimorphic_genes_pit_wide[sexually_dimorphic_genes_pit_wide$logFC > 0,]$gene
down_genes <- sexually_dimorphic_genes_pit_wide[sexually_dimorphic_genes_pit_wide$logFC < 0,]$gene


library(ggplot2)
library(dplyr)

# Convert top_genes_duncan to a data frame if it's not already
top_genes_alonso_df <- as.data.frame(top_genes_alonso)

# Add a column to indicate gene direction (up, down, or NS - not significant/not in our lists)
top_genes_alonso_df <- top_genes_alonso_df %>%
  mutate(gene_direction = case_when(
    rownames(.) %in% up_genes ~ "Up",
    rownames(.) %in% down_genes  ~ "Down",
    TRUE ~ "NS"
  ))



volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_alonso_df,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "alonso_volcano_plot_sex_de_min_3.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "alonso_volcano_plot_sex_de_min_3.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)




volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_alonso_df,show_ns= TRUE,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "alonso_volcano_plot_sex_de_min_3_ns.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "alonso_volcano_plot_sex_de_min_3_Ns.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)



##############
####Hou et al. (wt whole pituitaries) pituitary-wide validations!
####
##############

#First for genes that are sex-biased in at least one cell type

sexually_dimorphic_genes <- read.csv("/Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv")
sexually_dimorphic_genes_pit_wide <- sexually_dimorphic_genes[sexually_dimorphic_genes$occurs > 1,]
sexually_dimorphic_genes_pit_wide


up_genes <- sexually_dimorphic_genes_pit_wide[sexually_dimorphic_genes_pit_wide$logFC > 0,]$gene
down_genes <- sexually_dimorphic_genes_pit_wide[sexually_dimorphic_genes_pit_wide$logFC < 0,]$gene




# Convert top_genes_hou to a data frame
top_genes_hou_df <- as.data.frame(top_genes_hou)

# Assuming you want to color based on the same up_genes and down_genes lists
# Define gene direction based on presence in the lists and logFC direction
top_genes_hou_df <- top_genes_hou_df %>%
  mutate(gene_direction = case_when(
    rownames(.) %in% up_genes ~ "Up",
    rownames(.) %in% down_genes ~ "Down",
    TRUE ~ "NS"
  ))



volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_hou_df, y_gap_min = 9, y_gap_max = 10,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)



#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_1.png"), plot = volcano_plot_with_gap, width = 5, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_1.svg"), plot = volcano_plot_with_gap, width = 5, height = 3)


volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_hou_df, y_gap_min = 10, y_gap_max = 10,color1="#ffa500",color2="#63b3ed",highlight_gene = c("Fshb"))
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_1_fshb.png"), plot = volcano_plot_with_gap, width = 5, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_1_fshb.svg"), plot = volcano_plot_with_gap, width = 5, height = 3)




volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_hou_df, y_gap_min = 10, y_gap_max = 10, show_ns= TRUE,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_1_ns.png"), plot = volcano_plot_with_gap, width =5, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_1_ns.svg"), plot = volcano_plot_with_gap, width = 5, height = 3)





#First for genes that are sex-biased in at least three cell types


sexually_dimorphic_genes <- read.csv("/Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv")
sexually_dimorphic_genes_pit_wide <- sexually_dimorphic_genes[sexually_dimorphic_genes$occurs > 3,]
sexually_dimorphic_genes_pit_wide


up_genes <- sexually_dimorphic_genes_pit_wide[sexually_dimorphic_genes_pit_wide$logFC > 0,]$gene
down_genes <- sexually_dimorphic_genes_pit_wide[sexually_dimorphic_genes_pit_wide$logFC < 0,]$gene




# Convert top_genes_hou to a data frame
top_genes_hou_df <- as.data.frame(top_genes_hou)

# Assuming you want to color based on the same up_genes and down_genes lists
# Define gene direction based on presence in the lists and logFC direction
top_genes_hou_df <- top_genes_hou_df %>%
  mutate(gene_direction = case_when(
    rownames(.) %in% up_genes ~ "Up",
    rownames(.) %in% down_genes ~ "Down",
    TRUE ~ "NS"
  ))


volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_hou_df, y_gap_min = 10, y_gap_max = 10,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_3.png"), plot = volcano_plot_with_gap, width = 5, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_3.svg"), plot = volcano_plot_with_gap, width = 5, height = 3)




volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_hou_df, y_gap_min = 20, y_gap_max = 30, show_ns= TRUE,color1="#ffa500",color2="#63b3ed")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_3_ns.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_sex_de_min_3_ns.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)


