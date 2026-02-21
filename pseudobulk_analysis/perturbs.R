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

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)

#make boxplot of colsums log y
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")


exp_table = expression_table#[,meta_data$exptype == "sc"]
meta_data = meta_data#[meta_data$exptype == "sc",]


#prepare perturbations for the design matrix later


hfd = c("SRX21241800","SRX21241801", #Ruggiero-Ruff et al. (2024)
        "SRX21986106","SRX21986108","SRX21986113", "SRX21986114","SRX21986115","SRX21986116", #Miles et al. (2023)
        "SRX24151434","SRX24151435",#Qian et al. (2024)
        "SRX31166358","SRX31166359" #Guo et al. (2025)
        )


gata2_cko = c("SRX13290060",
              "SRX13290061",
              "SRX13290062")



#other perturbations: POMC, AX,GX, AX/GX, p25OE, Acvr1b KO, Tgfbr1 KO, LEPR KO, Stress, 





#### Assign variables
meta_data$prop1_df_df <- ifelse(meta_data$SRA_ID %in% prop1_df_df, "P", "C")
meta_data$tpit_status <- ifelse(meta_data$SRA_ID %in% tpit_ko, "KO",
                                ifelse(meta_data$SRA_ID %in% tpit_het, "HET", "WT"))
meta_data$hfd <- ifelse(meta_data$SRA_ID %in% hfd, "HFD", "C")

#for #
#for SRX24151438, SRX24151439 set as IREPKO
irepko = c("SRX24151438", "SRX24151439")

meta_data$hfd <- ifelse(meta_data$SRA_ID %in% irepko, "irepko", meta_data$hfd)

meta_data$hesx1_dn_rar <- ifelse(meta_data$SRA_ID %in% hesx1_dn_rar, "HESX1_DN_RAR", "C")
meta_data$gata2_cko <- ifelse(meta_data$SRA_ID %in% gata2_cko, "GATA2_CKO", "C")

#find what assignments there are where prop1_mut is M
print(unique(meta_data$assignments[meta_data$hfd  == "HFD"]))
#keep only these in meta_data and exp table
assignments_to_keep <- unique(meta_data$assignments[meta_data$hfd == "HFD"])

exp_table <- exp_table[,meta_data$assignments %in% assignments_to_keep]
meta_data <- meta_data[meta_data$assignments %in% assignments_to_keep,]

dge = DGEList(counts = exp_table, group = meta_data$assignments)

# Create the design matrix
design <- model.matrix(~0 + assignments + Author + assignments:hfd, data = meta_data)
#change to make.names
colnames(design) <- make.names(colnames(design))

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 500,min.prop= 0.2)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)

dge <- calcNormFactors(dge, method = "TMM")

library(variancePartition)
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~0 + assignments +  assignments:hfd   +  (1 | Author) +  (1 | Modality)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, meta_data, BPPARAM = param)

# fit dream model with contrasts
fit <- dream(vobjDream, form, meta_data)#, L)

head(coef(fit))

fit <- variancePartition::eBayes(fit, trend=FALSE)

head(coef(fit))

top_genes <- variancePartition::topTable(fit, coef="assignmentsStem_cells:hfdHFD", number = Inf, sort.by = "P",adjust.method="BH", lfc=0.25)
top_genes$genes <- rownames(top_genes)
print(head(top_genes,20))









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

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)

#make boxplot of colsums log y
boxplot(colSums(expression_table), main = "Boxplot of colSums of expression_table", ylab = "log10(colSums)", log = "y")


exp_table = expression_table#[,meta_data$exptype == "sc"]
meta_data = meta_data#[meta_data$exptype == "sc",]


#prepare perturbations for the design matrix later


hfd = c("SRX21241800","SRX21241801", #Ruggiero-Ruff et al. (2024)
        "SRX21986106","SRX21986108","SRX21986113", "SRX21986114","SRX21986115","SRX21986116", #Miles et al. (2023)
        "SRX24151434","SRX24151435",#Qian et al. (2024)
        "SRX31166358","SRX31166359" #Guo et al. (2025)
)


gata2_cko = c("SRX13290060",
              "SRX13290061",
              "SRX13290062")



meta_data$gata2_cko <- ifelse(meta_data$SRA_ID %in% gata2_cko, "GATA2_CKO", "C")

#find what assignments there are where prop1_mut is M
print(unique(meta_data$assignments[meta_data$gata2_cko  == "GATA2_CKO"]))
#keep only these in meta_data and exp table
assignments_to_keep <- unique(meta_data$assignments[meta_data$gata2_cko == "GATA2_CKO"])

exp_table <- exp_table[,meta_data$assignments %in% assignments_to_keep]
meta_data <- meta_data[meta_data$assignments %in% assignments_to_keep,]

dge = DGEList(counts = exp_table, group = meta_data$assignments)

# Create the design matrix
design <- model.matrix(~0 + assignments + Author + assignments:gata2_cko, data = meta_data)
#change to make.names
colnames(design) <- make.names(colnames(design))

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 500,min.prop= 0.2)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)

dge <- calcNormFactors(dge, method = "TMM")

library(variancePartition)
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~0 + assignments +  assignments:gata2_cko   +  (1 | Author) +  (1 | Modality)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, meta_data, BPPARAM = param)

# fit dream model with contrasts
fit <- dream(vobjDream, form, meta_data)#, L)

head(coef(fit))

fit <- variancePartition::eBayes(fit, trend=FALSE)

head(coef(fit))

top_genes <- variancePartition::topTable(fit, coef="assignmentsGonadotrophs:gata2_ckoGATA2_CKO", number = Inf, sort.by = "P",adjust.method="BH", lfc=0.25)
top_genes$genes <- rownames(top_genes)
print(head(top_genes,20))




