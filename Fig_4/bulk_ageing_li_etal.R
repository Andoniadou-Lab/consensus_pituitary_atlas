# Load required libraries
library(Matrix)
library(limma)
library(edgeR)
library(dplyr)

# Define sample paths and IDs
sample_dirs <- c(
  "SRX17837195",
  "SRX17837196",
  "SRX17837197",
  "SRX17837198",
  "SRX17837199",
  "SRX17837200",
  "SRX17837201",
  "SRX17837202"
)

base_path <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/bulk_ageing_li_et_al/samples/"

# Create a group factor (first 4 young, last 4 old)
group <- factor(c(rep("Young", 4), rep("Old", 4)))
print(data.frame(SampleID = sample_dirs, Group = group))

# Function to read a single-cell matrix from the files
sample_dir = "SRX17837202"
read_sc_matrix <- function(sample_dir) {
  full_path <- paste0(base_path, sample_dir, "/")
  
  # Read the Matrix Market file (.mtx)
  mat_file <- paste0(full_path, "cells_x_genes.total.mtx")
  mat <- readMM(mat_file)
  
  # Read gene names
  genes_file <- paste0(full_path, "cells_x_genes.genes.names.txt")
  genes <- readLines(genes_file)
  
  # Read barcodes (cell IDs)
  barcodes_file <- paste0(full_path, "cells_x_genes.barcodes.txt")
  barcodes <- readLines(barcodes_file)
  
  # Set dimension names
  colnames(mat) <- genes
  rownames(mat) <- barcodes
  
  #transposing the matrix
  mat <- t(mat)
  
  return(mat)
}

# Read matrices for each sample
sample_matrices <- list()
for (i in seq_along(sample_dirs)) {
  cat("Reading sample:", sample_dirs[i], "\n")
  sample_matrices[[i]] <- read_sc_matrix(sample_dirs[i])
}

# Combine counts by summing across cells for each sample
# This converts single-cell data to bulk RNA-seq-like data
combined_counts <- data.frame(row.names = rownames(sample_matrices[[1]]))

for (i in seq_along(sample_matrices)) {
  # Sum counts across all cells for each gene
  sample_sum <- rowSums(sample_matrices[[i]])
  combined_counts[[sample_dirs[i]]] <- sample_sum
}

# Check data structure
dim(combined_counts)
head(combined_counts)

# Create DGEList object
dge <- DGEList(counts = combined_counts, group = group)

# Filter out low-count genes
keep <- filterByExpr(dge, min.total.count=100,min.prop=0.5)
dge <- dge[keep, ]


dim(dge)

# Normalize
dge <- calcNormFactors(dge,method = "TMM")

# Design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
print(design)

# Voom transformation
v <- voom(dge, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design)

# Create contrast matrix for Young vs Old
cont.matrix <- makeContrasts(
  YoungVsOld = Old - Young,
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Get top DE genes
top_genes <- topTable(fit2, coef = "YoungVsOld", n = Inf, sort.by = "P")
top_genes <- top_genes %>% 
  arrange(adj.P.Val) %>%
  mutate(Gene = rownames(.)) %>%
  dplyr::select(Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# Save results
write.csv(top_genes, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/bulk_ageing_li_et_al/DE_results_young_vs_old.csv", row.names = FALSE)

# Show top 20 differentially expressed genes
top_20 <- top_genes[1:20, ]
print(top_20)

# Summary of significant genes
summary_sig <- summary(decideTests(fit2, p.value = 0.05))
print("Number of significantly DE genes (FDR < 0.05):")
print(summary_sig)

# Generate a volcano plot
with(top_genes, plot(logFC, -log10(P.Value), pch = 20, main = "Volcano Plot"))
with(subset(top_genes, adj.P.Val < 0.05), points(logFC, -log10(P.Value), pch = 20, col = "red"))
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)

# Generate an MD plot
limma::plotMD(fit2, column = 1, status = decideTests(fit2, p.value = 0.05)[,1], 
              main = "MD Plot: Young vs Old")





final_results = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/epitome/data/aging/v_0.01/aging_genes.csv")
#keep only padj < 0.05
final_results_pan_pit <- final_results[final_results$adj.P.Val < 0.05,]
#remove where cell_type "Mesenchymal_cells" "Pituicytes"   "Immune_cells" "Endothelial_cells"
final_results_pan_pit <- final_results_pan_pit[!final_results_pan_pit$cell_type %in% c("Pituicytes"), ]#"Mesenchymal_cells", "Pituicytes", "Immune_cells", "Endothelial_cells"),]
#calc avg logFC and geom avg adj pval
final_results_pan_pit$avg_logFC <- ave(final_results_pan_pit$logFC, final_results_pan_pit$genes, FUN = mean)
#geom avg adj pval
final_results_pan_pit$avg_adj_pval <- ave(final_results_pan_pit$adj.P.Val, final_results_pan_pit$genes, FUN = function(x) prod(x)^(1/length(x)))
#add another column saying number of times a given gene occurs
final_results_pan_pit$gene_count <- ave(final_results_pan_pit$genes, final_results_pan_pit$genes, FUN = length)
#keep only those where log2fc sign is consistent
final_results_pan_pit <- final_results_pan_pit[ave(final_results_pan_pit$logFC > 0, final_results_pan_pit$genes, FUN = function(x) length(unique(x))) == 1,]
final_results_pan_pit <- final_results_pan_pit[!duplicated(final_results_pan_pit$genes),]
#at least 3 gene count
final_results_pan_pit <- final_results_pan_pit[final_results_pan_pit$gene_count >= 3,]
final_results_pan_pit
#keep these genes final_results_pan_pit
genes_to_keep <- final_results_pan_pit$genes
#keep these in toptable
top_genes <- top_genes[top_genes$Gene %in% genes_to_keep,]
#recalc padj using bh
top_genes$adj.P.Val <- p.adjust(top_genes$P.Value, method = "BH")

#plot log2fc between top_genes and final_results_pan_pit
library(ggplot2)
merged <- merge(top_genes, final_results_pan_pit, by.x = "Gene", by.y = "genes")
ggplot(merged, aes(x = logFC.x, y = logFC.y)) + geom_point() 


#calc spearman
cor.test(merged$logFC.x, merged$avg_logFC, method = "spearman")
merged $color <- factor(ifelse(merged$logFC.x * merged$avg_logFC > 0, "gray", "blue"))

plot<- ggplot(merged, aes(x = logFC.x, y = avg_logFC, color = color)) +
  geom_point(aes(color = color), alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_color_manual(values = c("gray", "blue")) +  # Manually specify colors
  theme_minimal() + 
  theme(legend.position = "none") +
  xlab("Validation - log2FC") + 
  ylab("Discovery - avglog2FC")

figs_folder = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/bulk_ageing_li_et_al/figs/"
#make this folder if doesnt exist
if (!dir.exists(figs_folder)) {
  dir.create(figs_folder)
}

#save to figs folder /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/sox2_egfp_sorted/figs
ggsave(paste0(figs_folder, "ageing_comparison_scatter_logfc.png"), plot, width = 3, height = 3, dpi = 300)
#svg
ggsave(paste0(figs_folder, "ageing_comparison_scatter_logfc.svg"), plot, width = 3, height = 3, dpi = 300)

#what percentage is blue
sum(merged$color == "gray")/nrow(merged)



plot_consistent_volcano_with_gap <- function(df, y_gap_min = NULL, y_gap_max = NULL, show_ns = FALSE, p_val_threshold = 0.05,
                                             color1="blue",color2="purple"
) {
  
  consistent_sign_down <- rownames(df[(df$logFC > 0) &
                                        (df$gene_direction == "Down"),])
  consistent_sign_up <- rownames(df[(df$logFC < 0) &
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
  n_consistent_down_above_thresh <- sum((df$logFC > 0) & (df$gene_direction == "Down") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  total_down_above_thresh <- sum((df$gene_direction == "Down") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  
  n_consistent_up_above_thresh <- sum((df$logFC < 0) & (df$gene_direction == "Up") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  total_up_above_thresh <- sum((df$gene_direction == "Up") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  
  print(paste0("Number of consistent down-regulated genes above p-value threshold: ", n_consistent_down_above_thresh))
  print(paste0("Total down-regulated genes above p-value threshold: ", total_down_above_thresh))
  print(paste0("Number of consistent up-regulated genes above p-value threshold: ", n_consistent_up_above_thresh))
  print(paste0("Total up-regulated genes above p-value threshold: ", total_up_above_thresh))
  
  
  
  # Calculate max y value before trimming for annotation placement
  max_y <- max(-log10(df$adj.P.Val), na.rm = TRUE)
  
  # Create the initial volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    # Plot "NS" points first (in the background) if show_ns is TRUE
    {if (show_ns) geom_point(data = subset(df, gene_direction == "NS"),
                             aes(color = gene_direction), size = 1, alpha = 0.4)} +
    # Plot "Up" and "Down" points on top
    geom_point(data = subset(df, gene_direction %in% c("Up", "Down")),
               aes(color = gene_direction), size = 1, alpha = 0.4) +
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
             label = paste0(n_consistent_down, " (", n_consistent_down_above_thresh, ") consistent"), color = color1) +
    annotate("text", x = max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.2, hjust = 0.5,
             label = paste0("out of ", total_down, " (", total_down_above_thresh, ")"), color = color1) +
    
    # Annotate for Up-regulated genes
    annotate("text", x = -max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.3, hjust = 0.5,
             label = paste0(n_consistent_up, " (", n_consistent_up_above_thresh, ") consistent"), color = color2) +
    annotate("text", x = -max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.2, hjust = 0.5,
             label = paste0("out of ", total_up, " (", total_up_above_thresh, ")"), color =color2)
  
  # Add y-axis break if y_gap_min and y_gap_max are provided
  if (!is.null(y_gap_min) && !is.null(y_gap_max)) {
    p <- p + scale_y_break(c(y_gap_min, y_gap_max))
  }
  
  return(p)
}


top_genes_li_df <- as.data.frame(topTable(fit2, coef = "YoungVsOld", n = Inf, sort.by = "P")
                                 )
up_genes <- final_results_pan_pit[final_results_pan_pit$logFC > 0, "genes"]
down_genes <- final_results_pan_pit[final_results_pan_pit$logFC < 0, "genes"]

# Assuming you want to color based on the same up_genes and down_genes lists
# Define gene direction based on presence in the lists and logFC direction
top_genes_li_df <- top_genes_li_df %>%
  mutate(gene_direction = case_when(
    rownames(.) %in% up_genes ~ "Down",
    rownames(.) %in% down_genes ~ "Up",
    TRUE ~ "NS"
  ))

table(top_genes_li_df$gene_direction)



volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_li_df,color2 = "#FFA500", color1 = "#63B3ED")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_min3.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_min3.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)




volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_li_df,show_ns = TRUE,,color2 = "#FFA500", color1 = "#63B3ED")
print(volcano_plot_with_gap)

471+367
660+451
838/1111

(282+249) / (368+295)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_ns_min3.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_ns_min3.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)

(282+249)
(368+295)






#########
### Repeat for 2 months old comparison
########

final_results = read.csv( paste0("/Users/k23030440/epitome_code/epitome/data/aging/v_0.01/aging_genes_2_months.csv"))
#keep only padj < 0.05
final_results_pan_pit <- final_results[final_results$adj.P.Val < 0.05,]
#remove where cell_type "Mesenchymal_cells" "Pituicytes"   "Immune_cells" "Endothelial_cells"
final_results_pan_pit <- final_results_pan_pit[!final_results_pan_pit$cell_type %in% c("Pituicytes"), ]#"Mesenchymal_cells", "Pituicytes", "Immune_cells", "Endothelial_cells"),]
#calc avg logFC and geom avg adj pval
final_results_pan_pit$avg_logFC <- ave(final_results_pan_pit$logFC, final_results_pan_pit$genes, FUN = mean)
#geom avg adj pval
final_results_pan_pit$avg_adj_pval <- ave(final_results_pan_pit$adj.P.Val, final_results_pan_pit$genes, FUN = function(x) prod(x)^(1/length(x)))
#add another column saying number of times a given gene occurs
final_results_pan_pit$gene_count <- ave(final_results_pan_pit$genes, final_results_pan_pit$genes, FUN = length)
#keep only those where log2fc sign is consistent
final_results_pan_pit <- final_results_pan_pit[ave(final_results_pan_pit$logFC > 0, final_results_pan_pit$genes, FUN = function(x) length(unique(x))) == 1,]
final_results_pan_pit <- final_results_pan_pit[!duplicated(final_results_pan_pit$genes),]
#at least 3 gene count
final_results_pan_pit <- final_results_pan_pit[final_results_pan_pit$gene_count >= 3,]
final_results_pan_pit
#keep these genes final_results_pan_pit
genes_to_keep <- final_results_pan_pit$genes
#keep these in toptable
top_genes <- top_genes[top_genes$Gene %in% genes_to_keep,]
#recalc padj using bh
top_genes$adj.P.Val <- p.adjust(top_genes$P.Value, method = "BH")

#plot log2fc between top_genes and final_results_pan_pit
library(ggplot2)
merged <- merge(top_genes, final_results_pan_pit, by.x = "Gene", by.y = "genes")
ggplot(merged, aes(x = logFC.x, y = logFC.y)) + geom_point() 


#calc spearman
cor.test(merged$logFC.x, merged$avg_logFC, method = "spearman")
merged $color <- factor(ifelse(merged$logFC.x * merged$avg_logFC > 0, "gray", "blue"))

plot<- ggplot(merged, aes(x = logFC.x, y = avg_logFC, color = color)) +
  geom_point(aes(color = color), alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_color_manual(values = c("gray", "blue")) +  # Manually specify colors
  theme_minimal() + 
  theme(legend.position = "none") +
  xlab("Validation - log2FC") + 
  ylab("Discovery - avglog2FC")

figs_folder = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/bulk_ageing_li_et_al/figs/"
#make this folder if doesnt exist
if (!dir.exists(figs_folder)) {
  dir.create(figs_folder)
}

#save to figs folder /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/sox2_egfp_sorted/figs
ggsave(paste0(figs_folder, "ageing_comparison_scatter_logfc_2_months.png"), plot, width = 3, height = 3, dpi = 300)
#svg
ggsave(paste0(figs_folder, "ageing_comparison_scatter_logfc_2_months.svg"), plot, width = 3, height = 3, dpi = 300)

#what percentage is blue
sum(merged$color == "gray")/nrow(merged)



plot_consistent_volcano_with_gap <- function(df, y_gap_min = NULL, y_gap_max = NULL, show_ns = FALSE, p_val_threshold = 0.05,
                                             color1="blue",color2="purple"
) {
  
  consistent_sign_down <- rownames(df[(df$logFC > 0) &
                                        (df$gene_direction == "Down"),])
  consistent_sign_up <- rownames(df[(df$logFC < 0) &
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
  n_consistent_down_above_thresh <- sum((df$logFC > 0) & (df$gene_direction == "Down") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  total_down_above_thresh <- sum((df$gene_direction == "Down") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  
  n_consistent_up_above_thresh <- sum((df$logFC < 0) & (df$gene_direction == "Up") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  total_up_above_thresh <- sum((df$gene_direction == "Up") & (df$adj.P.Val < p_val_threshold), na.rm = TRUE)
  
  print(paste0("Number of consistent down-regulated genes above p-value threshold: ", n_consistent_down_above_thresh))
  print(paste0("Total down-regulated genes above p-value threshold: ", total_down_above_thresh))
  print(paste0("Number of consistent up-regulated genes above p-value threshold: ", n_consistent_up_above_thresh))
  print(paste0("Total up-regulated genes above p-value threshold: ", total_up_above_thresh))
  
  
  
  # Calculate max y value before trimming for annotation placement
  max_y <- max(-log10(df$adj.P.Val), na.rm = TRUE)
  
  # Create the initial volcano plot
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    # Plot "NS" points first (in the background) if show_ns is TRUE
    {if (show_ns) geom_point(data = subset(df, gene_direction == "NS"),
                             aes(color = gene_direction), size = 1, alpha = 0.4)} +
    # Plot "Up" and "Down" points on top
    geom_point(data = subset(df, gene_direction %in% c("Up", "Down")),
               aes(color = gene_direction), size = 1, alpha = 0.4) +
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
             label = paste0(n_consistent_down, " (", n_consistent_down_above_thresh, ") consistent"), color = color1) +
    annotate("text", x = max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.2, hjust = 0.5,
             label = paste0("out of ", total_down, " (", total_down_above_thresh, ")"), color = color1) +
    
    # Annotate for Up-regulated genes
    annotate("text", x = -max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.3, hjust = 0.5,
             label = paste0(n_consistent_up, " (", n_consistent_up_above_thresh, ") consistent"), color = color2) +
    annotate("text", x = -max(df$logFC, na.rm = TRUE) * 0.5, y = max_y * 0.2, hjust = 0.5,
             label = paste0("out of ", total_up, " (", total_up_above_thresh, ")"), color =color2)
  
  # Add y-axis break if y_gap_min and y_gap_max are provided
  if (!is.null(y_gap_min) && !is.null(y_gap_max)) {
    p <- p + scale_y_break(c(y_gap_min, y_gap_max))
  }
  
  return(p)
}


top_genes_li_df <- as.data.frame(topTable(fit2, coef = "YoungVsOld", n = Inf, sort.by = "P")
)
up_genes <- final_results_pan_pit[final_results_pan_pit$logFC > 0, "genes"]
down_genes <- final_results_pan_pit[final_results_pan_pit$logFC < 0, "genes"]

# Assuming you want to color based on the same up_genes and down_genes lists
# Define gene direction based on presence in the lists and logFC direction
top_genes_li_df <- top_genes_li_df %>%
  mutate(gene_direction = case_when(
    rownames(.) %in% up_genes ~ "Down",
    rownames(.) %in% down_genes ~ "Up",
    TRUE ~ "NS"
  ))

table(top_genes_li_df$gene_direction)



volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_li_df,color2 = "#FFA500", color1 = "#63B3ED")
print(volcano_plot_with_gap)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_min3_2_months.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_min3_2_months.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)




volcano_plot_with_gap <- plot_consistent_volcano_with_gap(top_genes_li_df,show_ns = TRUE,,color2 = "#FFA500", color1 = "#63B3ED")
print(volcano_plot_with_gap)

471+367
660+451
838/1111

(282+249) / (368+295)

#save plot to figs_folder
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_ns_min3.png"), plot = volcano_plot_with_gap, width = 4, height = 3)
#svg
ggsave(paste0(figs_folder, "volcano_plot_age_de_stem_li_ns_min3.svg"), plot = volcano_plot_with_gap, width = 4, height = 3)

(282+249)
(368+295)










