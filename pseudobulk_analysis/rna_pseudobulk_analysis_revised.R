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
prop1_index = which(rownames(expression_table)=="Prop1")
prop1_sum <- sum(expression_table[prop1_index,])
keep_genes <- filterByExpr(dge, design, min.total.count = 500,min.prop= 0.5)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)
keep_genes["Stat4"]
dge_sc <- dge_sc[gene_names,]
dge_sn <- dge_sn[gene_names,]
dge_multi <- dge_multi[gene_names,]
dge_parse <- dge_parse[gene_names,]

#print genes names with top counts
rownames(dge)[order(rowSums(dge$counts), decreasing = TRUE)[1:15]]
#also print their total count
rowSums(dge$counts)[order(rowSums(dge$counts), decreasing = TRUE)[1:15]]


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

fit <- voom(dge, design, plot=TRUE)


meta_data$assay_version <- meta_data$`10X version`

form <- ~ (1 | Modality) + (1 | Sex) + (1 | Author)
vp <- fitExtractVarPartModel(fit, form, meta_data)

plotVarPart(sortCols(vp))




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
      levels = c("Modality", "Author", "Sex", "Residuals")
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
ggsave(paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/revision_figures/rna_variance_partition_boxplot_",version,".png"), width = 3.5, height = 2.5)
#svg
ggsave(paste0("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/figures/revision_figures/rna_variance_partition_boxplot_",version,".svg"), width = 3.5, height = 2.5)






#double running this based on
#https://support.bioconductor.org/p/114663/
#first changing the weights to the separate ones


corfit <- duplicateCorrelation(fit, design, block=meta_data$exptype)
corfit$consensus.correlation
fit <- voom(dge, design, plot=TRUE, block=meta_data$exptype, correlation=corfit$consensus.correlation)
corfit <- duplicateCorrelation(fit, design, block=meta_data$exptype)
corfit$consensus.correlation

fit <- lmFit(fit, design, block=meta_data$exptype, correlation=corfit$consensus.correlation)





#____________
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~0 + assignments + (1 | Author) + (1 | Modality)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, meta_data, BPPARAM = param)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation

fit <- dream(vobjDream, form, meta_data)
fit <- eBayes(fit)


###
# Make column names valid
colnames(design) <- make.names(colnames(design), unique = TRUE)

print(head(coef(fit)))

#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
write.csv(as.data.frame(coef(fit)), "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/coef_2026_01_20.csv")
coefs_table <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/coef_2026_01_20.csv")
#save to epitome

write.csv(as.data.frame(coef(fit)), "/Users/k23030440/epitome_code/epitome/data/heatmap/v_0.02/coef.csv")

#turn coef(fit) into a dataframe
coef_df = as.data.frame(coef(fit))

coef_df$gene = rownames(coef_df)

# Create contrast - first lets try for Stem cell comparisons

L <- makeContrastsDream(form, meta_data,
                        contrasts = c( 
                          
                          #Corticotrophs
                          assignmentsCorticotrophs__assignmentsGonadotrophs = "assignmentsCorticotrophs - assignmentsGonadotrophs",
                          assignmentsCorticotrophs__assignmentsMelanotrophs = "assignmentsCorticotrophs - assignmentsMelanotrophs",
                          assignmentsCorticotrophs__assignmentsThyrotrophs = "assignmentsCorticotrophs - assignmentsThyrotrophs",
                          assignmentsCorticotrophs__assignmentsSomatotrophs = "assignmentsCorticotrophs - assignmentsSomatotrophs",
                          assignmentsCorticotrophs__assignmentsLactotrophs = "assignmentsCorticotrophs - assignmentsLactotrophs",
                          assignmentsCorticotrophs__assignmentsPituicytes = "assignmentsCorticotrophs - assignmentsPituicytes",
                          assignmentsCorticotrophs__assignmentsImmune_cells = "assignmentsCorticotrophs - assignmentsImmune_cells",
                          assignmentsCorticotrophs__assignmentsMesenchymal_cells = "assignmentsCorticotrophs - assignmentsMesenchymal_cells",
                          assignmentsCorticotrophs__assignmentsEndothelial_cells = "assignmentsCorticotrophs - assignmentsEndothelial_cells",
                          assignmentsCorticotrophs__assignmentsStem_cells = "assignmentsCorticotrophs - assignmentsStem_cells",
                          
                          #Gonadotrophs
                          assignmentsGonadotrophs__assignmentsCorticotrophs = "assignmentsGonadotrophs - assignmentsCorticotrophs",
                          assignmentsGonadotrophs__assignmentsMelanotrophs = "assignmentsGonadotrophs - assignmentsMelanotrophs",
                          assignmentsGonadotrophs__assignmentsThyrotrophs = "assignmentsGonadotrophs - assignmentsThyrotrophs",
                          assignmentsGonadotrophs__assignmentsSomatotrophs = "assignmentsGonadotrophs - assignmentsSomatotrophs",
                          assignmentsGonadotrophs__assignmentsLactotrophs = "assignmentsGonadotrophs - assignmentsLactotrophs",
                          assignmentsGonadotrophs__assignmentsPituicytes = "assignmentsGonadotrophs - assignmentsPituicytes",
                          assignmentsGonadotrophs__assignmentsImmune_cells = "assignmentsGonadotrophs - assignmentsImmune_cells",
                          assignmentsGonadotrophs__assignmentsMesenchymal_cells = "assignmentsGonadotrophs - assignmentsMesenchymal_cells",
                          assignmentsGonadotrophs__assignmentsEndothelial_cells = "assignmentsGonadotrophs - assignmentsEndothelial_cells",
                          assignmentsGonadotrophs__assignmentsStem_cells = "assignmentsGonadotrophs - assignmentsStem_cells",
                          
                          #Melanotrophs
                          assignmentsMelanotrophs__assignmentsGonadotrophs = "assignmentsMelanotrophs - assignmentsGonadotrophs",
                          assignmentsMelanotrophs__assignmentsCorticotrophs = "assignmentsMelanotrophs - assignmentsCorticotrophs",
                          assignmentsMelanotrophs__assignmentsThyrotrophs = "assignmentsMelanotrophs - assignmentsThyrotrophs",
                          assignmentsMelanotrophs__assignmentsSomatotrophs = "assignmentsMelanotrophs - assignmentsSomatotrophs",
                          assignmentsMelanotrophs__assignmentsLactotrophs = "assignmentsMelanotrophs - assignmentsLactotrophs",
                          assignmentsMelanotrophs__assignmentsPituicytes = "assignmentsMelanotrophs - assignmentsPituicytes",
                          assignmentsMelanotrophs__assignmentsImmune_cells = "assignmentsMelanotrophs - assignmentsImmune_cells",
                          assignmentsMelanotrophs__assignmentsMesenchymal_cells = "assignmentsMelanotrophs - assignmentsMesenchymal_cells",
                          assignmentsMelanotrophs__assignmentsEndothelial_cells = "assignmentsMelanotrophs - assignmentsEndothelial_cells",
                          assignmentsMelanotrophs__assignmentsStem_cells = "assignmentsMelanotrophs - assignmentsStem_cells",
                          
                          
                          #Thyrotrophs
                          assignmentsThyrotrophs__assignmentsGonadotrophs = "assignmentsThyrotrophs - assignmentsGonadotrophs",
                          assignmentsThyrotrophs__assignmentsCorticotrophs = "assignmentsThyrotrophs - assignmentsCorticotrophs",
                          assignmentsThyrotrophs__assignmentsMelanotrophs = "assignmentsThyrotrophs - assignmentsMelanotrophs",
                          assignmentsThyrotrophs__assignmentsSomatotrophs = "assignmentsThyrotrophs - assignmentsSomatotrophs",
                          assignmentsThyrotrophs__assignmentsLactotrophs = "assignmentsThyrotrophs - assignmentsLactotrophs",
                          assignmentsThyrotrophs__assignmentsPituicytes = "assignmentsThyrotrophs - assignmentsPituicytes",
                          assignmentsThyrotrophs__assignmentsImmune_cells = "assignmentsThyrotrophs - assignmentsImmune_cells",
                          assignmentsThyrotrophs__assignmentsMesenchymal_cells = "assignmentsThyrotrophs - assignmentsMesenchymal_cells",
                          assignmentsThyrotrophs__assignmentsEndothelial_cells = "assignmentsThyrotrophs - assignmentsEndothelial_cells",
                          assignmentsThyrotrophs__assignmentsStem_cells = "assignmentsThyrotrophs - assignmentsStem_cells",
                          
                          
                          
                          #Somatotrophs
                          assignmentsSomatotrophs__assignmentsGonadotrophs = "assignmentsSomatotrophs - assignmentsGonadotrophs",
                          assignmentsSomatotrophs__assignmentsCorticotrophs = "assignmentsSomatotrophs - assignmentsCorticotrophs",
                          assignmentsSomatotrophs__assignmentsMelanotrophs = "assignmentsSomatotrophs - assignmentsMelanotrophs",
                          assignmentsSomatotrophs__assignmentsThyrotrophs = "assignmentsSomatotrophs - assignmentsThyrotrophs",
                          assignmentsSomatotrophs__assignmentsLactotrophs = "assignmentsSomatotrophs - assignmentsLactotrophs",
                          assignmentsSomatotrophs__assignmentsPituicytes = "assignmentsSomatotrophs - assignmentsPituicytes",
                          assignmentsSomatotrophs__assignmentsImmune_cells = "assignmentsSomatotrophs - assignmentsImmune_cells",
                          assignmentsSomatotrophs__assignmentsMesenchymal_cells = "assignmentsSomatotrophs - assignmentsMesenchymal_cells",
                          assignmentsSomatotrophs__assignmentsEndothelial_cells = "assignmentsSomatotrophs - assignmentsEndothelial_cells",
                          assignmentsSomatotrophs__assignmentsStem_cells = "assignmentsSomatotrophs - assignmentsStem_cells",
                          
                          
                          #Lactotrophs
                          assignmentsLactotrophs__assignmentsGonadotrophs = "assignmentsLactotrophs - assignmentsGonadotrophs",
                          assignmentsLactotrophs__assignmentsCorticotrophs = "assignmentsLactotrophs - assignmentsCorticotrophs",
                          assignmentsLactotrophs__assignmentsMelanotrophs = "assignmentsLactotrophs - assignmentsMelanotrophs",
                          assignmentsLactotrophs__assignmentsThyrotrophs = "assignmentsLactotrophs - assignmentsThyrotrophs",
                          assignmentsLactotrophs__assignmentsSomatotrophs = "assignmentsLactotrophs - assignmentsSomatotrophs",
                          assignmentsLactotrophs__assignmentsPituicytes = "assignmentsLactotrophs - assignmentsPituicytes",
                          assignmentsLactotrophs__assignmentsImmune_cells = "assignmentsLactotrophs - assignmentsImmune_cells",
                          assignmentsLactotrophs__assignmentsMesenchymal_cells = "assignmentsLactotrophs - assignmentsMesenchymal_cells",
                          assignmentsLactotrophs__assignmentsEndothelial_cells = "assignmentsLactotrophs - assignmentsEndothelial_cells",
                          assignmentsLactotrophs__assignmentsStem_cells = "assignmentsLactotrophs - assignmentsStem_cells",
                          
                          
                          #Pituicytes
                          assignmentsPituicytes__assignmentsGonadotrophs = "assignmentsPituicytes - assignmentsGonadotrophs",
                          assignmentsPituicytes__assignmentsCorticotrophs = "assignmentsPituicytes - assignmentsCorticotrophs",
                          assignmentsPituicytes__assignmentsMelanotrophs = "assignmentsPituicytes - assignmentsMelanotrophs",
                          assignmentsPituicytes__assignmentsThyrotrophs = "assignmentsPituicytes - assignmentsThyrotrophs",
                          assignmentsPituicytes__assignmentsSomatotrophs = "assignmentsPituicytes - assignmentsSomatotrophs",
                          assignmentsPituicytes__assignmentsLactotrophs = "assignmentsPituicytes - assignmentsLactotrophs",
                          assignmentsPituicytes__assignmentsImmune_cells = "assignmentsPituicytes - assignmentsImmune_cells",
                          assignmentsPituicytes__assignmentsMesenchymal_cells = "assignmentsPituicytes - assignmentsMesenchymal_cells",
                          assignmentsPituicytes__assignmentsEndothelial_cells = "assignmentsPituicytes - assignmentsEndothelial_cells",
                          assignmentsPituicytes__assignmentsStem_cells = "assignmentsPituicytes - assignmentsStem_cells",
                          
                          
                          #Immune_cells
                          assignmentsImmune_cells__assignmentsGonadotrophs = "assignmentsImmune_cells - assignmentsGonadotrophs",
                          assignmentsImmune_cells__assignmentsCorticotrophs = "assignmentsImmune_cells - assignmentsCorticotrophs",
                          assignmentsImmune_cells__assignmentsMelanotrophs = "assignmentsImmune_cells - assignmentsMelanotrophs",
                          assignmentsImmune_cells__assignmentsThyrotrophs = "assignmentsImmune_cells - assignmentsThyrotrophs",
                          assignmentsImmune_cells__assignmentsSomatotrophs = "assignmentsImmune_cells - assignmentsSomatotrophs",
                          assignmentsImmune_cells__assignmentsLactotrophs = "assignmentsImmune_cells - assignmentsLactotrophs",
                          assignmentsImmune_cells__assignmentsPituicytes = "assignmentsImmune_cells - assignmentsPituicytes",
                          assignmentsImmune_cells__assignmentsMesenchymal_cells = "assignmentsImmune_cells - assignmentsMesenchymal_cells",
                          assignmentsImmune_cells__assignmentsEndothelial_cells = "assignmentsImmune_cells - assignmentsEndothelial_cells",
                          assignmentsImmune_cells__assignmentsStem_cells = "assignmentsImmune_cells - assignmentsStem_cells",
                          
                          
                          #Mesenchymal_cells
                          assignmentsMesenchymal_cells__assignmentsGonadotrophs = "assignmentsMesenchymal_cells - assignmentsGonadotrophs",
                          assignmentsMesenchymal_cells__assignmentsCorticotrophs = "assignmentsMesenchymal_cells - assignmentsCorticotrophs",
                          assignmentsMesenchymal_cells__assignmentsMelanotrophs = "assignmentsMesenchymal_cells - assignmentsMelanotrophs",
                          assignmentsMesenchymal_cells__assignmentsThyrotrophs = "assignmentsMesenchymal_cells - assignmentsThyrotrophs",
                          assignmentsMesenchymal_cells__assignmentsSomatotrophs = "assignmentsMesenchymal_cells - assignmentsSomatotrophs",
                          assignmentsMesenchymal_cells__assignmentsLactotrophs = "assignmentsMesenchymal_cells - assignmentsLactotrophs",
                          assignmentsMesenchymal_cells__assignmentsPituicytes = "assignmentsMesenchymal_cells - assignmentsPituicytes",
                          assignmentsMesenchymal_cells__assignmentsImmune_cells = "assignmentsMesenchymal_cells - assignmentsImmune_cells",
                          assignmentsMesenchymal_cells__assignmentsEndothelial_cells = "assignmentsMesenchymal_cells - assignmentsEndothelial_cells",
                          assignmentsMesenchymal_cells__assignmentsStem_cells = "assignmentsMesenchymal_cells - assignmentsStem_cells",
                          
                          
                          #Endothelial_cells
                          assignmentsEndothelial_cells__assignmentsGonadotrophs = "assignmentsEndothelial_cells - assignmentsGonadotrophs",
                          assignmentsEndothelial_cells__assignmentsCorticotrophs = "assignmentsEndothelial_cells - assignmentsCorticotrophs",
                          assignmentsEndothelial_cells__assignmentsMelanotrophs = "assignmentsEndothelial_cells - assignmentsMelanotrophs",
                          assignmentsEndothelial_cells__assignmentsThyrotrophs = "assignmentsEndothelial_cells - assignmentsThyrotrophs",
                          assignmentsEndothelial_cells__assignmentsSomatotrophs = "assignmentsEndothelial_cells - assignmentsSomatotrophs",
                          assignmentsEndothelial_cells__assignmentsLactotrophs = "assignmentsEndothelial_cells - assignmentsLactotrophs",
                          assignmentsEndothelial_cells__assignmentsPituicytes = "assignmentsEndothelial_cells - assignmentsPituicytes",
                          assignmentsEndothelial_cells__assignmentsImmune_cells = "assignmentsEndothelial_cells - assignmentsImmune_cells",
                          assignmentsEndothelial_cells__assignmentsMesenchymal_cells = "assignmentsEndothelial_cells - assignmentsMesenchymal_cells",
                          assignmentsEndothelial_cells__assignmentsStem_cells = "assignmentsEndothelial_cells - assignmentsStem_cells",
                          
                          
                          #Stem cells
                          assignmentsStem_cells__assignmentsGonadotrophs = "assignmentsStem_cells - assignmentsGonadotrophs",
                          assignmentsStem_cells__assignmentsCorticotrophs = "assignmentsStem_cells - assignmentsCorticotrophs",
                          assignmentsStem_cells__assignmentsMelanotrophs = "assignmentsStem_cells - assignmentsMelanotrophs",
                          assignmentsStem_cells__assignmentsThyrotrophs = "assignmentsStem_cells - assignmentsThyrotrophs",
                          assignmentsStem_cells__assignmentsSomatotrophs = "assignmentsStem_cells - assignmentsSomatotrophs",
                          assignmentsStem_cells__assignmentsLactotrophs = "assignmentsStem_cells - assignmentsLactotrophs",
                          assignmentsStem_cells__assignmentsPituicytes = "assignmentsStem_cells - assignmentsPituicytes",
                          assignmentsStem_cells__assignmentsImmune_cells = "assignmentsStem_cells - assignmentsImmune_cells",
                          assignmentsStem_cells__assignmentsMesenchymal_cells = "assignmentsStem_cells - assignmentsMesenchymal_cells",
                          assignmentsStem_cells__assignmentsEndothelial_cells = "assignmentsStem_cells - assignmentsEndothelial_cells"
 
                        )
)


# fit dream model with contrasts
fit <- dream(vobjDream, form, meta_data, L)
fit <- eBayes(fit)

top_genes <- topTable(fit, coef="assignmentsGonadotrophs__assignmentsThyrotrophs", number = Inf, sort.by = "P",adjust.method="BH", lfc=0.25)
top_genes$genes <- rownames(top_genes)
print(head(top_genes,20))


top_genes <- top_genes %>% dplyr::mutate(abs_logFC = abs(logFC))

#hist
hist(top_genes$abs_logFC, breaks=50, main="lfc distribution", xlab="abs logFC")

#how many sig
top_genes %>% filter(adj.P.Val < 0.05)

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
    

    
    contrast_name_final <- paste0(contrast_name, "__", other_contrast_name)

    
    sig_genes <- topTable(fit, coef=contrast_name_final, number = Inf, sort.by = "P",adjust.method="BH", lfc=0.5) %>%
      rownames_to_column("gene") %>%
      filter(adj.P.Val < 0.05 & logFC > log2fc) %>%
      dplyr::mutate(group1 = celltype,
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
write.csv(all_comparisons, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/all_comparisons_2026_01_21.csv")

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
  print(n_comps)
  print(groupings1)
  print(groupings2)
  
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
    
    print(sig_genes)
    
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
  
  #remove duplicate rows
  sig_genes_list <- sig_genes_list %>% distinct()
  
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
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% dplyr::mutate(mean_log2fc = mean(log2fc))
  #add geom_mean_adj_pval
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% dplyr::mutate(geom_mean_adj_pval = exp(mean(log(pvalue))))
  #add mean AvgExpr
  sig_genes_list <- sig_genes_list %>% group_by(gene) %>% dplyr::mutate(mean_AvgExpr = mean(AveExpr))
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

#/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/epitome/data/gene_group_annotation/v_0.02/cpdb.csv
cpdb <- read_csv("/Users/k23030440/epitome_code/epitome/data/gene_group_annotation/v_0.01/cpdb.csv")
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
all_results <- all_results %>% dplyr::mutate(TF = ifelse(gene %in% tf_list, 1, 0))
table(all_results$direction,all_results$grouping)

#save to /Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/
write_csv(all_results, "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/grouping_lineage_markers_0909.csv")
write_csv(all_results, paste0("/Users/k23030440/epitome_code/epitome/data/markers/",version,"/grouping_lineage_markers.csv"))
write_csv(all_results, paste0("/Users/k23030440/epitome_code/epitome/data/heatmap/",version,"/rna_grouping_lineage_markers.csv"))

#load lineage markers
all_results <- read_csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/grouping_lineage_markers_0909.csv")
table(all_results$grouping,all_results$direction)

library(edgeR)
library(Matrix)

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
expression_table_original <- as(expression_table_original, "CsparseMatrix")

normalized_expression <- normalize_expression(expression_table_original, dge_original)

normalized_expression <- as(normalized_expression, "CsparseMatrix")

#first 5x5 values
normalized_expression[1:5,1:5]

version = "v_0.01"
Matrix::writeMM(normalized_expression, paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/normalized_data.mtx"))
write_csv(meta_data_original, paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/meta_data.csv"))
write.table(rownames(normalized_expression), paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(normalized_expression), paste0("/Users/k23030440/epitome_code/epitome/data/expression/",version,"/datasets.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

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

#remove multi_rna
#expression_table <- expression_table[,meta_data$exptype != "multi_rna"]
#meta_data <- meta_data[meta_data$exptype != "multi_rna",]



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
meta_data$sample = colnames(expression_table)


expression_table <- expression_table[,order(meta_data$exptype)]
meta_data <- meta_data[order(meta_data$exptype),]


#print all values
print(meta_data$exptype)
#add a new column called index
meta_data$index <- colnames(expression_table)
#make index the rownames
rownames(meta_data) <- meta_data$index
expression_table[c(1:5),c(1:5)]
meta_data$psbulk_n_cells[c(1:5)]
expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")
expression_table[c(1:5),c(1:5)]

expression_table <- expression_table * median(meta_data$psbulk_n_cells)
expression_table[c(1:5),c(1:5)]





exp_table_sc = expression_table[,meta_data$exptype == "sc"]
meta_data_sc = meta_data[meta_data$exptype == "sc",]

exp_table_sn = expression_table[,meta_data$exptype == "sn"]
meta_data_sn = meta_data[meta_data$exptype == "sn",]


#keep only gonadotrophs
gonadotrophs_exp_table = expression_table[,meta_data$assignments == "Gonadotrophs"]
gonadotrophs_meta_data = meta_data[meta_data$assignments == "Gonadotrophs",]

#plot Gal for male and female
gal_index = which(rownames(gonadotrophs_exp_table) == "Gal")
male_index = which(gonadotrophs_meta_data$Comp_sex==1)
female_index = which(gonadotrophs_meta_data$Comp_sex==0)

male_gal = gonadotrophs_exp_table[gal_index,male_index]
male_gal
gonadotrophs_meta_data$SRA_ID[male_index]

female_gal = gonadotrophs_exp_table[gal_index,female_index]
female_gal
gonadotrophs_meta_data$SRA_ID[female_index]

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
#plot(ecdf(rowSums(dge$counts)), main = "ecdf of total counts of genes", xlab = "Total counts", ylab = "Proportion of genes")
#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 500,min.prop= 0.2)
#look at Gpr101
keep_genes[rownames(dge) == "Megf11"]
print(sum(keep_genes))
dge <- dge[keep_genes, ]
dim(dge)
gene_names <- rownames(dge)
dge_sc <- dge_sc[gene_names,]
dge_sn <- dge_sn[gene_names,]
dge_multi <- dge_multi[gene_names,]
#dge_parse <- dge_parse[gene_names,]



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

colnames(design)<-make.names(colnames(design))

#double running this based on
#https://support.bioconductor.org/p/114663/


fit <- voom(dge, design, plot=TRUE)
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
diff_expression <- eBayes(diff_expression)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf,lfc=1, adjust.method="BH")
top_genes$genes <- rownames(top_genes)
print(rownames((top_genes)))
#how many sig
print(nrow(top_genes %>% filter(adj.P.Val < 0.05 )))

cpdb <- read_csv("/Users/k23030440/epitome_code/epitome/data/gene_group_annotation/v_0.02//cpdb.csv")
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
  results <- topTable(diff_expr, number=Inf,adjust.method="BH", lfc=1)
  
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
  dplyr::mutate(total_diff_exp = sum(abs(count))) %>% # Total absolute differentially expressed genes
  ungroup() %>%
  arrange(desc(total_diff_exp)) %>%            # Order by total genes
  dplyr::mutate(cell_type = factor(cell_type, levels = unique(cell_type))) # Maintain order

# Add a column for fixed label positions
plot_data_long <- plot_data_long %>%
  dplyr::mutate(label_position = ifelse(count > 0, 400, -400)) # Labels at fixed positions

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
# Save the plot as PNG and SVG
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
  results <- topTable(diff_expr, number=Inf,adjust.method="BH")
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
#only save those with abs log2fc > 1
df_final <- df %>% filter(abs(logFC) > 1)
df_final  <- df_final  %>%
  select(-occurs)

df_final <- df_final  %>%
  group_by(gene) %>%
  summarize(occurs = n_distinct(cell_type)) %>%
  left_join(df_final , by = "gene")

df_final 
write.csv(df_final,"/Users/k23030440/epitome_code/epitome/data/sex_dimorphism/v_0.02/sexually_dimorphic_genes.csv")

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
#keep where logfc > 1
dim(df_sex_genes_filtered )
df_sex_genes_filtered <- df_sex_genes_filtered %>% filter(abs(logFC) > 1)
df_sex_genes_filtered <- df_sex_genes_filtered %>%
  select(-occurs)
df_sex_genes_filtered  <- df_sex_genes_filtered   %>%
  group_by(gene) %>%
  summarize(occurs = n_distinct(cell_type)) %>%
  left_join(df_sex_genes_filtered, by = "gene")

df_sex_genes_filtered <-df_sex_genes_filtered %>% filter(occurs > 4)

#keep only those where sign is same direction in all cases
df_sex_genes_filtered <- df_sex_genes_filtered %>%
  group_by(gene) %>%
  filter(all(logFC > 0) | all(logFC < 0)) %>%
  ungroup()

#how many unique genes
print(df_sex_genes_filtered %>% distinct(gene) %>% nrow())


"Gal" %in%  df_sex_genes_filtered$gene
df_sex_genes <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes_raw.csv")
#keep genes in df_sex_filtered
df_sex_genes <- df_sex_genes %>% filter(gene %in% df_sex_genes_filtered$gene)

df_sex_genes <- df_sex_genes %>%
  group_by(gene) %>%
  dplyr::mutate(avg_log2fc = mean(logFC, na.rm = TRUE)) %>%
  ungroup()

#first 20 gene names
genes_of_interest <- unique(df_sex_genes$gene)

final_results <- df_sex_genes
"Gal" %in% genes_of_interest

# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)


#only keep those that are avg abs logfc > 1
final_results <- final_results %>%
  group_by(gene) %>%
  dplyr::mutate(avg_abs_log2fc = mean(abs(logFC), na.rm = TRUE)) %>%
  ungroup()

#remove genes starting with gm or ens !grepl("^Gm|^Gm[0-9]+$|^ENS" !grepl("Rik$"
final_results <- final_results %>% filter(!grepl("^Gm|^ENS|Rik$", gene))

#keep top 30 genes with highest abs logFC
final_results <- final_results %>% arrange(desc(avg_abs_log2fc))
genes_of_interest1 <- unique(final_results$gene)[1:30]

final_results<- final_results %>% filter(gene %in% genes_of_interest1)

final_results <- final_results %>% arrange(desc(avg_log2fc))
genes_of_interest <- unique(final_results$gene)[1:30]

#remove fshb
genes_of_interest <- genes_of_interest[genes_of_interest != "Fshb"]






#-1 the logFC
final_results$logFC <- -final_results$logFC
# Create the heatmap
plot <- final_results %>%
  # Filter to only include genes of interest
  dplyr::filter(gene %in% genes_of_interest) %>%
  # Create a significance flag and preserve input order
  dplyr::mutate(is_significant = adj.P.Val < 0.05,
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
                       linewidth = 2
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









