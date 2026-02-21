library(DESeq2)
library(Seurat)
library(tidyverse)
library(reticulate)
library(variancePartition)


version="v_0.02"
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

# Create the Seurat object
seurat_object <- CreateSeuratObject(counts = counts, meta.data = meta_data)
#print unique Author
rownames(seurat_object) <- vars
seurat_object$assignments <- seurat_object$new_cell_type

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

#remove parse
expression_table <- expression_table[,meta_data$exptype != "parse"]
meta_data <- meta_data[meta_data$exptype != "parse",]

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

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

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
keep_genes <- filterByExpr(dge, design, min.total.count = 500,min.prop= 0.5)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)
dge_sc <- dge_sc[gene_names,]
dge_sn <- dge_sn[gene_names,]
dge_multi <- dge_multi[gene_names,]
dge_parse <- dge_parse[gene_names,]


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

fit <- voom(dge, design, plot=TRUE)

form <- ~  (1 | exptype) + (1 | Sex)  + (1 | Author) 
vp <- fitExtractVarPartModel(fit, form, meta_data)

plotVarPart(sortCols(vp))





######------------------------------------------
# Repeat for ATAC
######------------------------------------------




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

form <- ~  (1 | exptype) + (1 | Sex) + (1 | Age_numeric) +(1 | Author) 
vp <- fitExtractVarPartModel(fit, form, meta_data)

plotVarPart(sortCols(vp))


#def function

plot_variance_partition <- function(vp, version) {
  #Boxplot with ggplot
  vp_long <- as.data.frame(vp)
  
  library(tidyverse)
  
  # Pivot to long format
  vp_plot <- vp_long %>%
    rownames_to_column("feature") %>%  # keep genomic coordinate if needed
    pivot_longer(
      cols = -feature,
      names_to = "Coefficient",
      values_to = "Value"
    )
  
  vp_plot <- vp_plot %>%
    mutate(
      Coefficient = factor(
        Coefficient,
        levels = c("Age_numeric", "Author", "exptype", "Sex", "Residuals")
      )
    )
  
  # Nature-style boxplot
  ggplot(vp_plot, aes(x = Coefficient, y = Value)) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,     # NO outlier dots
      fill = "white",
      color = "black",
      linewidth = 0.4
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11, color = "black"),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) + theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )) +
    labs(
      x = NULL,
      y = "Variance explained"
    )
  #save to /Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/figures/revision_figures
  ggsave(paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/revision_figures/atac_variance_partition_boxplot_",version,".png"), width = 3.5, height = 2.5)
  #svg
  ggsave(paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/revision_figures/atac_variance_partition_boxplot_",version,".svg"), width = 3.5, height = 2.5)
}

