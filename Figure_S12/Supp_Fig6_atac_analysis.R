if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("annotatr"),force=TRUE)

library(reticulate)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(AnnotationDbi)
library(annotatr)
library(IRanges) # Often loaded with GenomicRanges, but good to explicitly load
library(GenomicRanges) # Core package for GRanges objects

mem.maxVSize(1e+10)


#install
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", force=TRUE)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#loading libraries
library(DESeq2)
library(limma)
library(edgeR)
library(tidyverse)
library(org.Mm.eg.db)
library(Seurat)




#keep only those rows that are in
all_results <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/grouping_lineage_markers.csv")
#keep obly sig
all_results <- all_results[all_results$geom_mean_adj_pval<0.05,]
#keep where direction is up and grouping is grouping_1
all_results <- all_results[all_results$direction == "up" & all_results$grouping == "grouping_1",]

peaks_to_keep = all_results$gene






######Repeating by extracting actula motif positions

library(Seurat)
library(Signac)

# Specify the path to the .h5ad file
#h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pb_h5ad.h5ad"
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

#keep those in motifs$grouping == grouping_1_up
#seurat_object <- seurat_object[all_results$gene, ]

root = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/"
figs_folder = paste(root,"Figures/",sep="")


meta_data = seurat_object@meta.data
meta_data$index <- 1:nrow(meta_data)
#rename sample to GEO
meta_data$GEO = meta_data$sample
print(length(unique(meta_data$Author)))
print(length(unique(meta_data$SRA_ID)))



#load meta2 from /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pituitary_atlas_atac.xlsx
#meta2 = readxl::read_excel("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pituitary_atlas_atac.xlsx")
#load meta2 /Users/k23030440/Downloads/pituitary_atlas_atac.xlsx
meta2 = readxl::read_excel("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/pituitary_atlas (3).xlsx")
#add missing columns to meta_data from meta2. shared columns is GEO. keep the same row dimensions
#iterate through meta_data
# First make meta2 have only first occurrence of each GEO
meta2_first <- meta2[!duplicated(meta2$GEO), ]

# Then merge will maintain row dimensions
meta_data <- merge(meta_data, meta2_first, by="GEO", all.x=TRUE)
#reorder based on index
meta_data <- meta_data[order(meta_data$index),]

seurat_object <- AddMetaData(seurat_object, meta_data)

##############
#Signac enrichment
##############
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(Signac)
matrix = seurat_object[["RNA"]]$counts
gc()
#make it matrix
matrix <- as(matrix, "sparseMatrix")

gc()
atac_assay <- CreateChromatinAssay(matrix,sep = c(":", "-"),)

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





motif_names <- unname(unlist(atac_object@assays$ATAC@motifs@motif.names))

#find motifs that are relevant /Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/tfs_retrieve_from_epitome/results/multimodal_hits_combined.csv
motifs = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tfs_retrieve_from_epitome/results/multimodal_hits_combined.csv")
#keep only where motifs motifs$grouping == grouping_1
motifs = motifs[motifs$grouping == "grouping_1_up",]
motif_genes = motifs$gene
#caps them
motif_genes = toupper(motif_genes)

#upper all motif_names
motif_names = toupper(motif_names)

#which entries in motif_names contain on of motifs$gene
indices = c()
genes = c()
for (i in 1:length(motif_genes)) {
  #find where it equals motifnames
  indices_i = which(motif_names == motif_genes[i])
  indices = c(indices, indices_i)
  genes = c(genes, motif_genes[i])
  print(unname(indices_i))
  print(motif_genes[i])
}

motifs_sub_matrix = atac_object@assays$ATAC@motifs@data[, indices]


#keep only those rows that are in
all_results <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/grouping_lineage_markers.csv")
#keep obly sig
all_results <- all_results[all_results$geom_mean_adj_pval<0.05,]
#keep where direction is up and grouping is grouping_1
all_results <- all_results[all_results$direction == "up" & all_results$grouping == "grouping_1",]
#in all_results$gene change : to -
all_results$gene <- gsub(":", "-", all_results$gene)

motif_names = atac_object@assays$ATAC@motifs@motif.names
#find motifs that are relevant /Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/tfs_retrieve_from_epitome/results/multimodal_hits_combined.csv
motifs = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/tfs_retrieve_from_epitome/results/multimodal_hits_combined.csv")
#keep only where motifs motifs$grouping == grouping_1
motifs = motifs[motifs$grouping == "grouping_1_up",]
motif_genes = motifs$gene
#caps them
motif_genes = toupper(motif_genes)
#caps motif_names
motif_names = toupper(motif_names)
#which entries in motif_names contain on of motifs$gene
indices = c()
genes = c()
for (i in 1:length(motif_genes)) {
  #find where it equals motifnames
  indices_i = which(motif_names == motif_genes[i])
  indices = c(indices, indices_i)
  genes = c(genes, motif_genes[i])
  print(unname(indices_i))
  print(motif_genes[i])
}

annotations_res <- c()
motifs_res <- c()
motif = c()

library(reticulate)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(AnnotationDbi)
library(annotatr)
library(IRanges) # Often loaded with GenomicRanges, but good to explicitly load
library(GenomicRanges) # C


total_occurrences_df = data.frame(Motif=character(), Total_Occurrences=integer(), stringsAsFactors=FALSE)

# Build the annotations
annots <- c('mm10_basicgenes','mm10_genes_intergenic', 'mm10_cpgs',  'mm10_genes_cds', 'mm10_genes_firstexons', 'mm10_genes_intronexonboundaries', 'mm10_genes_exonintronboundaries', 'mm10_lncrna_gencode', 'mm10_enhancers_fantom')

annotations.mm10 <- build_annotations(genome = 'mm10', annotations = annots)

for (i in 1:length(indices)) {
  
  index = indices[i]
  
  
  
  # Get the column name
  
  # Get the column name
  motif_name <- colnames(motifs_sub_matrix)[i]
  print(motif_name)
  
  
  #change the first - to : but keep the second -
  peaks_to_keep <- gsub("^(chr)?(\\w+)-(\\d+)-(\\d+)$", "\\1\\2:\\3-\\4", peaks_to_keep)
  
  sig_peaks<- GRanges(seqnames = sub("^(chr)?", "chr", sub(":(\\d+)-(\\d+)", "", peaks_to_keep)), # Ensure chromosome names start with "chr"
                     ranges = IRanges(start = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\1", peaks_to_keep)), 
                                      end = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\2", peaks_to_keep))),
                     strand = "*")
  
  motif_name <- colnames(atac_object@assays$ATAC@motifs@data)[index]
  #motif_name = "MA0077.1"
  positions <- atac_object@assays$ATAC@motifs@positions[motif_name]
  single_motif_granges <- positions[[motif_name]]
  #concat seqnames and ranges to a peak format
  peaks<- paste0(
    as.character(seqnames(single_motif_granges)),
    ":",
    start(single_motif_granges),
    "-",
    end(single_motif_granges)
  )
  
  #change the first - to : but keep the second -
  peaks <- gsub("^(chr)?(\\w+)-(\\d+)-(\\d+)$", "\\1\\2:\\3-\\4", peaks)
  
  my_peaks<- GRanges(seqnames = sub("^(chr)?", "chr", sub(":(\\d+)-(\\d+)", "", peaks)), # Ensure chromosome names start with "chr"
                     ranges = IRanges(start = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\1", peaks)), 
                                      end = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\2", peaks))),
                     strand = "*") # Assuming no strand information is available
  
  my_peaks <- subsetByOverlaps(x = my_peaks, ranges = sig_peaks)
  print(length(my_peaks))
  print("Actually keeping")
  my_peaks <- subsetByOverlaps(x = my_peaks, ranges = annotations.mm10)
  print(length(my_peaks))
  
  seqlevelsStyle(my_peaks) <- "UCSC" # Ensure chromosome names match UCSC style (e.g., "chr1")
  seqinfo(my_peaks) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(my_peaks)]
  
  #remove duplicate peaks
  my_peaks <- unique(my_peaks)
  
  # Annotate your peaks
  annotated_peaks <- annotate_regions(
    regions = my_peaks,
    annotations = annotations.mm10,
    minoverlap = 3,
    ignore.strand = TRUE
  )
  
  annotated_peaks
  
  
  
  
  #keep first instance of each peak
  #annotated_peaks <- annotated_peaks[!duplicated(ranges(annotated_peaks$annot)), ]
  
  
  annotation_summary <- summarize_annotations(
    annotated_peaks
  )
  annotation_summary
  #add another row with annot.type total
  annotation_summary <- rbind(annotation_summary, data.frame(annot.type="total", n=length(my_peaks)))
  annotation_summary["motif"]<- genes[i]
  annotation_summary
  annotations_res[[length(annotations_res) + 1]] <- annotation_summary
  motifs_res <- c(motifs_res, genes[i])
  total_occurrences_df <- rbind(total_occurrences_df, data.frame(Motif=motif_genes[i], Total_Occurrences=length(my_peaks), stringsAsFactors=FALSE))
  motif <- c(motif, motif_genes[i])
}

# Combine all valid annotation summaries into a single data frame
all_annotation_summaries <- do.call(rbind, annotations_res)
all_annotation_summaries

# Use 'n' and 'annot.type' as confirmed from your output
annotation_matrix_df <- aggregate(n ~ motif + annot.type, data = all_annotation_summaries, sum)

# Cast to wide format (matrix)
library(tidyr)
annotation_matrix <- pivot_wider(annotation_matrix_df,
                                 names_from = annot.type, # Use annot.type here
                                 values_from = n,        # Use n here
                                 values_fill = 0)

# Convert to a matrix and set row names
motif_names_in_matrix <- annotation_matrix$motif
annotation_matrix <- as.matrix(annotation_matrix[,-1])
rownames(annotation_matrix) <- motif_names_in_matrix

# Handle cases where all values in a column might be zero (can cause issues with dist)
annotation_matrix <- annotation_matrix[, colSums(annotation_matrix) > 0, drop = FALSE]

# Handle cases where there might be only one motif left after filtering
if (nrow(annotation_matrix) < 2) {
  stop("Not enough motifs with valid annotations to perform hierarchical clustering. At least 2 motifs are required.")
}
if (ncol(annotation_matrix) < 1) {
  stop("No common annotation categories found for the selected motifs. Cannot create heatmap.")
}

# Create the heatmap with dendrogram
#normalise by total peaks
row_normalized_matrix <- annotation_matrix / annotation_matrix[,"total"]

#remove column total
row_normalized_matrix <- row_normalized_matrix[ , !(colnames(row_normalized_matrix) %in% c("total"))]

#log2
row_normalized_matrix <- log2(row_normalized_matrix + 1e-4) # Adding a small constant to avoid log(0)
#now z score for each col
row_normalized_matrix <- apply(row_normalized_matrix, 2, function(x) (x - mean(x)) / sd(x))

# Perform hierarchical clustering on the motifs (rows)
motif_dist <- dist(row_normalized_matrix, method = "euclidean")

motif_cluster <- hclust(motif_dist, method = "ward.D2")


#remove mm10_genes_introns
#row_normalized_matrix_final <- row_normalized_matrix[,- which(colnames(row_normalized_matrix) == "mm10_genes_introns"), drop = FALSE]

# Create the heatmap with the row-normalized data
#viridis cmap
library(viridis)
library(circlize)
library(ComplexHeatmap)
bwr_col <- colorRamp2(
  c(-4, 0, 4), # Breaks: min, midpoint (0), max
  c("blue", "white", "red") # Colors: Blue, White, Red
)

# 2. Update the Heatmap call
h <- Heatmap(row_normalized_matrix,
             name = "Annotation Count",
             cluster_rows = motif_cluster,
             cluster_columns = FALSE,
             show_row_dend = TRUE,
             show_column_dend = FALSE,
             # === CHANGED: Use the BWR color function ===
             col = bwr_col, 
             # ==========================================
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             heatmap_legend_param = list(title = "Z-score") # NOTE: Changed title to reflect Z-scores
)

h


#save /Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/figures/atac_fold_enrichment.svg and png with dpi 500
library(svglite)

# Save as SVG
svglite::svglite("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/atac_fold_enrichment.svg", width = 4, height = 6)
draw(h)
dev.off()

# Save as PNG with DPI of 500
png("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/atac_fold_enrichment.png", width = 4, height = 6, units = "in", res = 500)
draw(h)
dev.off()




# Assuming 'atac_object', 'indices', 'motif_genes', and 'BSgenome.Mmusculus.UCSC.mm10'
# are already defined in your environment from the previous steps.

library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)

# Get the TSS locations
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)

# Create an empty list to store TSS distances for each motif
tss_distances_by_motif <- list()

for (i in 1:length(indices)) {
  index = indices[i]
  # Get the column name
  motif_name <- colnames(atac_object@assays$ATAC@motifs@data)[index]
  positions <- atac_object@assays$ATAC@motifs@positions[motif_name]
  single_motif_granges <- positions[[motif_name]]
  
  print(length(single_motif_granges))
  single_motif_granges <- subsetByOverlaps(x = single_motif_granges, ranges = sig_peaks)
  print(length(single_motif_granges))
  # Ensure motif GRanges have UCSC style seqlevels for TSS distance calculation
  seqlevelsStyle(single_motif_granges) <- "UCSC"
  seqinfo(single_motif_granges) <- seqinfo(genes)[seqlevels(single_motif_granges)]
  
  # Calculate TSS distances
  # Find the nearest gene to each motif peak
  nearest_genes <- suppressWarnings(distanceToNearest(single_motif_granges, genes))
  
  # Get the actual distances
  distances <- mcols(nearest_genes)$distance
  
  # Store the distances
  tss_distances_by_motif[[motif_genes[i]]] <- distances
  
}

# Prepare data for plotting
plot_data <- do.call(rbind, lapply(names(tss_distances_by_motif), function(motif_name) {
  data.frame(
    Motif = motif_name,
    TSS_Distance = tss_distances_by_motif[[motif_name]]
  )
}))

# Calculate the percentage of values less than 10000 for each motif
percentage_data <- plot_data %>%
  group_by(Motif) %>%
  summarise(
    pct_less_than_10000 = mean(TSS_Distance < 5000) * 100
  ) %>%
  # Arrange motifs by percentage in descending order
  arrange(desc(pct_less_than_10000))

# Set the order of Motif levels in plot_data based on the arranged percentage_data
plot_data$Motif <- factor(plot_data$Motif, levels = percentage_data$Motif)


y_min_stagger <- 15
y_max_stagger <- 25

# Create a sequence of y positions for the labels
# This creates an alternating pattern
stagger_levels_original <- c(y_min_stagger, y_max_stagger)
num_motifs <- nrow(percentage_data)
percentage_data$staggered_y_position <- stagger_levels_original[((1:num_motifs - 1) %% length(stagger_levels_original)) + 1]


#add to tss distance a noise of mean 5 with sd 1
set.seed(42) # For reproducibility
plot_data$TSS_Distance <- plot_data$TSS_Distance + rnorm(nrow(plot_data), mean = 5, sd = 1)

# Create the stripplot
p <- ggplot(plot_data, aes(x = Motif, y = TSS_Distance)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.05, size = 0.2,
             #color "#0000ff"
             color = "#0000ff"
             ) +
  labs(
    title = "TSS Distances for Motif-Associated Peaks",
    x = "TF (ordered by % of TSS distance < 5kb)",
    y = "Distance to Nearest TSS (bp)"
  ) +
  theme_minimal() +
  # make y log10
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Add text labels for percentages
  geom_text(
    data = percentage_data, # Use the ordered percentage_data with staggered y
    aes(x = Motif, y = staggered_y_position, label = paste0(round(pct_less_than_10000, 0), "%")),
    vjust = 0, # Align text bottom to the staggered_y_position
    hjust = 0.5, # Center horizontally
    size = 2.7, # Increased text size for better readability
    inherit.aes = FALSE # Do not inherit aesthetics from the main plot_data
  ) +
  #add line with #0000ff at y = 10000
  geom_hline(yintercept = 5000, linetype = "dashed", color = "#0000ff")

p

#save
ggsave(
  filename = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/tss_distance_motif_peaks.png",
  plot = p,
  width = 8,
  height = 3.5,
  dpi = 600
)

#svg
svg(
  filename = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/tss_distance_motif_peaks.svg",
  width = 8,
  height = 3.5
)







































#now on all peaks 
my_peaks = rownames(atac_assay)
#change first - to : but keep the second -
my_peaks <- gsub("^(chr)?(\\w+)-(\\d+)-(\\d+)$", "\\1\\2:\\3-\\4", my_peaks)
my_peaks<- GRanges(seqnames = sub("^(chr)?", "chr", sub(":(\\d+)-(\\d+)", "", my_peaks)), # Ensure chromosome names start with "chr"
                   ranges = IRanges(start = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\1", my_peaks)), 
                                    end = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\2", my_peaks))),
                   strand = "*") # Assuming no strand information is available

seqlevelsStyle(my_peaks) <- "UCSC" # Ensure chromosome names match UCSC style (e.g., "chr1")
seqinfo(my_peaks) <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(my_peaks)]

annots <- c('mm10_basicgenes','mm10_genes_intergenic')

# Build the annotations
annotations.mm10 <- build_annotations(genome = 'mm10', annotations = annots)

# Annotate your peaks
annotated_peaks <- annotate_regions(
  regions = my_peaks,
  annotations = annotations.mm10,
  minoverlap = 1L,
  ignore.strand = TRUE
)


annotation_summary <- summarize_annotations(
  annotated_peaks
)

# Print the summary
print(annotation_summary)














# Load necessary libraries
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(memes)
library(universalmotif)
library(dplyr) # For %>% and arrange/slice

# 1. Read and filter the data as provided by the user
all_results <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE_atac/grouping_lineage_markers.csv")


all_results <- all_results[all_results$geom_mean_adj_pval < 0.05,]
table(all_results$direction, all_results$grouping)

all_results <- all_results[all_results$direction == "up" & all_results$grouping == "grouping_1",]
all_results$gene <- gsub(":", "-", all_results$gene)
#take top 1000 rows


# 2. Convert the 'gene' column (peak coordinates) into a GRanges object
# The 'gene' column is in "chr-start-end" format.
# We need to parse it to create a GRanges object.
peak_ranges_df <- data.frame(
  chrom = sapply(strsplit(all_results$gene, "-"), `[`, 1),
  start = as.integer(sapply(strsplit(all_results$gene, "-"), `[`, 2)),
  end = as.integer(sapply(strsplit(all_results$gene, "-"), `[`, 3))
)

# Create GRanges object
my_granges <- GRanges(
  seqnames = peak_ranges_df$chrom,
  ranges = IRanges(start = peak_ranges_df$start, end = peak_ranges_df$end)
)

# Resize GRanges to a fixed width (e.g., 100bp) centered on the original peak
# This is a common practice for motif discovery to focus on a central region.
#my_granges_resized <- GenomicRanges::resize(my_granges, width = 100, fix = "center")

# 3. Extract DNA sequences using the specified genome
genome <- BSgenome.Mmusculus.UCSC.mm10
peak_sequences <- get_sequence(my_granges, genome = genome)
names(peak_sequences) <- paste0("peak", seq_along(peak_sequences))

# 4. Perform de novo motif discovery using STREME
# Using shuffled input sequences as background.
# 4. Perform de novo motif discovery using STREME
# Using shuffled input sequences as background.
streme_results <- runStreme(
  input = peak_sequences,
  control = "shuffle",
  outdir = "streme_output_mm10_peaks", # Output directory
  alph = "dna",
  nmotifs = 10, # Discover top 3 motifs
  minw = 6,
  maxw = 30
)

# 5. Identify and display the top THREE enriched motifs
cat("--- Top 3 Enriched de Novo Motifs ---\n")

top_three_motifs_info <- streme_results %>%
  dplyr::arrange(pval) %>%
  dplyr::slice(1:3) # Get the top three rows

top_three_motifs <- as.list(top_three_motifs_info$motif) # Ensure it's a plain list

for (i in 1:length(top_three_motifs)) {
  current_motif <- top_three_motifs[[i]]
  current_motif_info_row <- top_three_motifs_info[i, ]
  
  cat(paste0("\nMotif ", i, ": ", current_motif_info_row$name, "\n"))
  cat("  Consensus: ", current_motif@consensus, "\n")
  cat("  P-value: ", current_motif_info_row$pval, "\n")
  cat("  E-value: ", current_motif_info_row$evalue, "\n")
  
  # Visualize each motif
  universalmotif::view_motifs(current_motif)
}




# 6. Infer motif identity by cross-checking against JASPAR (Known Motifs Database)
cat("\n--- Inferring Motif Identity using JASPAR ---\n")

jaspar_motifs_mm10 <- TFBSTools::getMatrixSet(
  x=JASPAR2022,
  opts = list(
    tax_group = "vertebrates",
    collection = "CORE",
    all_versions = FALSE
  )
)

#convert jaspar_motfs_to universalmotif
jaspar = convert_motifs(jaspar_motifs_mm10)
# Run Tomtom using the path to the temporary MEME file
tomtom_results <- runTomTom(
  input = streme_results,
  database = jaspar, # Target database
  outdir = "tomtom_output_mm10_peaks",
  thresh = 0.05
)

tomtom_results 

# You can then inspect 'tomtom_results' to see the matching motifs.
cat("\n--- TomTom Results (Top Matches) ---\n")
if (!is.null(tomtom_results) && nrow(tomtom_results) > 0) {
  top_matches <- tomtom_results %>%
    dplyr::group_by(query_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  print(top_matches)
} else {
  cat("No significant TomTom matches found.\n")
}

# (Optional) Clean up the temporary file when you're done
# file.remove(temp_meme_file)
