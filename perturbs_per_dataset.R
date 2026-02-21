## Perturbation effects
library(DESeq2)
library(Seurat)
library(tidyverse)
#update R

library(reticulate)
# Specify the path to the .h5ad file
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/pdatas_2026_02_14.h5ad"

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


gata2_cko = c("SRX13290060",
              "SRX13290061",
              "SRX13290062")

gata2_controls = c("SRX13290057",
                   "SRX13290058",
                   "SRX13290059")


#only keep first cko and control
expression_table <- expression_table[,meta_data$SRA_ID %in% c(gata2_cko, gata2_controls)]
meta_data <- meta_data[meta_data$SRA_ID %in% c(gata2_cko, gata2_controls),]

meta_data$gata2_cko <- ifelse(meta_data$SRA_ID %in% gata2_cko, "GATA2_CKO", "C")


#find what assignments there are where prop1_mut is M
assignments_that_occur_twice = table(meta_data$assignments)
assignments_that_occur_twice = names(assignments_that_occur_twice[assignments_that_occur_twice > 4])
assignments_to_keep <-assignments_that_occur_twice


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

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)

exp_table = expression_table
meta_data = meta_data

exp_table <- exp_table[,meta_data$assignments %in% assignments_to_keep]
meta_data <- meta_data[meta_data$assignments %in% assignments_to_keep,]

dge = DGEList(counts = exp_table, group = meta_data$assignments)

# Create the design matrix
design <- model.matrix(~0 + assignments + assignments:gata2_cko, data = meta_data)
#change to make.names
colnames(design) <- make.names(colnames(design))

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 50,min.prop= 0.2)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)

dge <- calcNormFactors(dge, method = "TMM")

vobj <- voom(dge, design, plot=TRUE)
fit <- lmFit(vobj, design)

fit <- eBayes(fit)

head(coef(fit))

top_genes <- topTable(fit, coef="assignmentsGonadotrophs.gata2_ckoGATA2_CKO", number = Inf, sort.by = "P",adjust.method="BH", lfc=0.5)
top_genes$genes <- rownames(top_genes)
print(head(top_genes,20))

coef_names <- c("assignmentsGonadotrophs.gata2_ckoGATA2_CKO", 
                "assignmentsMelanotrophs.gata2_ckoGATA2_CKO",
                "assignmentsSomatotrophs.gata2_ckoGATA2_CKO",
                "assignmentsLactotrophs.gata2_ckoGATA2_CKO",
                "assignmentsStem_cells.gata2_ckoGATA2_CKO")


#save perturb results to csv
perturbs_folder = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/"
for (coef in coef_names){
  print(coef)
  top_genes <- topTable(fit, coef=coef, number = Inf, sort.by = "P",adjust.method="BH", lfc=0.5)
  top_genes$genes <- rownames(top_genes)
  print(head(top_genes,20))
  write.csv(top_genes, paste(perturbs_folder,coef,"_GATA2_CKO_vs_Control_limma3v3.csv",sep=""))
}






library(tidyverse)

perturbs_folder <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/"

                
                
coef_names <- c(
  "assignmentsGonadotrophs.gata2_ckoGATA2_CKO", 
                "assignmentsMelanotrophs.gata2_ckoGATA2_CKO",
                "assignmentsSomatotrophs.gata2_ckoGATA2_CKO",
                "assignmentsLactotrophs.gata2_ckoGATA2_CKO",
                "assignmentsStem_cells.gata2_ckoGATA2_CKO")

all_results_gata2<- map_df(coef_names, function(coef){
  
  file <- paste0(perturbs_folder, coef, "_GATA2_CKO_vs_Control_limma3v3.csv")
  
  read.csv(file) %>%
    mutate(cell_type = coef)
})


#remove genes with Gm
#all_results_ablation <- all_results_ablation %>% filter(!str_detect(genes, "^Gm"))


sig_counts <- all_results_gata2 %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(cell_type = str_remove(cell_type, "assignments"),
         cell_type = str_remove(cell_type, ".gata2_ckoGATA2_CKO")) %>%
  count(cell_type, name = "n_sig")

# include zero-hit cell types
all_cell_types <- tibble(cell_type = coef_names) %>%
  mutate(cell_type = str_remove(cell_type, "assignments"),
         cell_type = str_remove(cell_type, ".gata2_ckoGATA2_CKO"))

sig_counts <- all_cell_types %>%
  left_join(sig_counts, by = "cell_type") %>%
  mutate(n_sig = replace_na(n_sig, 0))

p_ablation <- ggplot(sig_counts,
                     aes(x = reorder(cell_type, n_sig),
                         y = n_sig)) +
  geom_bar(stat = "identity", fill = "#0000ff", width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(x = NULL,
       y = "Number of significant DE genes",
       title = "DT ablation-induced differential expression") +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p_ablation

ggsave(paste0(perturbs_folder, "gata2_DE_gene_counts.png"),
       p_ablation, width = 4, height = 3, dpi = 300)

ggsave(paste0(perturbs_folder, "gata2_DE_gene_counts.svg"),
       p_ablation, width = 4, height = 3)













####
###------------------------------------------------------------------------------
# Cell ablation
####



library(DESeq2)
library(Seurat)
library(tidyverse)
#update R
library(reticulate)
# Specify the path to the .h5ad file
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/pdatas_2026_02_04.h5ad"

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

#only keep entries where Author is in
author_list = c("Vennekens et al. (2021)",
                "Laporte et al. (2022)",
                "De Vriendt et al. (2025)")

expression_table <- expression_table[,meta_data$Author %in% author_list]
meta_data <- meta_data[meta_data$Author %in% author_list,]


dt = c("ERX4978538","ERX4978541","ERX4978542","ERX4978543",#Vennekens et al. (2021) 9-11 weeks and 14 months
       "ERX7593482","ERX7593483", #Laporte et al. (2022) #7 days
       "ERX13202613","ERX13202614"#De Vriendt et al. (2025) #24-26 month
)


meta_data$dt <- ifelse(meta_data$SRA_ID %in% dt, "DT", "C")

#merge study and age
meta_data$study_age <- paste(meta_data$Author, meta_data$Age, sep="_")
#turn to make names
meta_data$study_age <- make.names(meta_data$study_age)

#where dt=dt and study_age is "Laporte.et.al...2022._7.days", change to young_dt
meta_data$dt <- ifelse(meta_data$dt == "DT" & meta_data$study_age == "Laporte.et.al...2022._7.days", "young_dt", meta_data$dt)

#find what assignments there are where prop1_mut is M
print(unique(meta_data$assignments[meta_data$dt  == "DT"]))
#keep only these in meta_data and exp table
assignments_to_keep <- unique(meta_data$assignments[meta_data$dt == "DT"])

table(meta_data$dt)

#find what assignments there are where prop1_mut is M
assignments_that_occur_twice = table(meta_data$assignments)
assignments_that_occur_twice = names(assignments_that_occur_twice[assignments_that_occur_twice > 6])
assignments_to_keep <-assignments_that_occur_twice


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

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)

exp_table = expression_table
meta_data = meta_data

exp_table <- exp_table[,meta_data$assignments %in% assignments_to_keep]
meta_data <- meta_data[meta_data$assignments %in% assignments_to_keep,]

dge = DGEList(counts = exp_table, group = meta_data$assignments)

# Create the design matrix



meta_data$dt_age <- case_when(
  meta_data$dt == "C" ~ "C",
  meta_data$dt == "DT" ~ "old_DT",      # From Vennekens & De Vriendt
  meta_data$dt == "young_dt" ~ "young_DT" # From Laporte
)
meta_data$dt_age <- factor(meta_data$dt_age, levels = c("C", "old_DT", "young_DT"))

# Pool the "old" studies together
meta_data$study_pooled <- ifelse(
  meta_data$Author == "Laporte et al. (2022)", 
  "Laporte", 
  "Vennekens_DeVriendt"
)
meta_data$study_pooled <- factor(meta_data$study_pooled)

# Check distribution
table(meta_data$assignments, meta_data$dt_age)
table(meta_data$study_pooled, meta_data$dt_age)



design <- model.matrix(~0 + assignments + assignments:dt_age, data = meta_data)
colnames(design) <- make.names(colnames(design))

# Check rank
cat("Rank:", qr(design)$rank, "Columns:", ncol(design), "\n")

design <- model.matrix(~0 + assignments + study_pooled + dt_age + assignments:dt_age, data = meta_data)
colnames(design) <- make.names(colnames(design))

# Check rank
cat("Rank:", qr(design)$rank, "Columns:", ncol(design), "\n")

#if DT or young_dt set to DT
meta_data$final_dt = ifelse(meta_data$dt == "DT" | meta_data$dt == "young_dt", "DT", "C")
meta_data$final_dt

design <- model.matrix(~0 + assignments:study_pooled + assignments:final_dt, data = meta_data)
#change to make.names
colnames(design) <- make.names(colnames(design))

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 50,min.prop= 0.1)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)

dge <- calcNormFactors(dge, method = "TMM")


vobj <- voom(dge, design, plot=TRUE)
fit <- lmFit(vobj, design)

fit <- eBayes(fit)

head(coef(fit))

top_genes <- topTable(fit, coef="assignmentsStem_cells.final_dtDT", number = Inf, sort.by = "P",adjust.method="BH")
top_genes$genes <- rownames(top_genes)
print(head(top_genes,20))

#contrast young and old DT effect
cont_matrix <- makeContrasts( final_dtDT - final_dtC,
  levels = design
)

fit2 <- contrasts.fit(fit, cont_matrix)
diff_expression <- eBayes(fit2)

# Get top differentially expressed genes
top_genes <- topTable(diff_expression, number = Inf, adjust.method="BH")
top_genes$genes <- rownames(top_genes)
top_genes 
print(rownames((top_genes)))


library(variancePartition)
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ 0 +  assignments + assignments:final_dt + (1 | Author)
print(dim(meta_data))
print(dim(exp_table))
library(stats)
# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, meta_data, BPPARAM = param)

# fit dream model with contrasts
fit <- dream(vobjDream, form, meta_data)#, L)

head(coef(fit))

fit <- variancePartition::eBayes(fit)

head(coef(fit))

top_genes <- topTable(fit, coef="assignmentsStem_cells:final_dtDT", number = Inf, sort.by = "P",adjust.method="BH")
top_genes$genes <- rownames(top_genes)
print(head(top_genes,20))



coef_names <- c("assignmentsCorticotrophs:final_dtDT", "assignmentsEndothelial_cells:final_dtDT", "assignmentsGonadotrophs:final_dtDT",
                
                "assignmentsLactotrophs:final_dtDT", "assignmentsMesenchymal_cells:final_dtDT", "assignmentsSomatotrophs:final_dtDT",
                
                "assignmentsStem_cells:final_dtDT", "assignmentsThyrotrophs:final_dtDT")


#save perturb results to csv
perturbs_folder = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/"
for (coef in coef_names){
  print(coef)
  top_genes <- topTable(fit, coef=coef, number = Inf, sort.by = "P",adjust.method="BH", lfc=0.25)
  top_genes$genes <- rownames(top_genes)
  print(head(top_genes,20))
  write.csv(top_genes, paste(perturbs_folder,coef,"_ablation_dream.csv",sep=""))
}


library(tidyverse)

perturbs_folder <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/"

coef_names <- c(
  "assignmentsCorticotrophs:final_dtDT",
  "assignmentsEndothelial_cells:final_dtDT",
  "assignmentsGonadotrophs:final_dtDT",
  "assignmentsLactotrophs:final_dtDT",
  "assignmentsMesenchymal_cells:final_dtDT",
  "assignmentsSomatotrophs:final_dtDT",
  "assignmentsStem_cells:final_dtDT",
  "assignmentsThyrotrophs:final_dtDT"
)

all_results_ablation <- map_df(coef_names, function(coef){
  
  file <- paste0(perturbs_folder, coef, "_ablation_dream.csv")
  
  read.csv(file) %>%
    mutate(cell_type = coef)
})


#remove genes with Gm
all_results_ablation <- all_results_ablation %>% filter(!str_detect(genes, "^Gm"))



sig_counts <- all_results_ablation %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(cell_type = str_remove(cell_type, "assignments"),
         cell_type = str_remove(cell_type, ":final_dtDT")) %>%
  count(cell_type, name = "n_sig")

# include zero-hit cell types
all_cell_types <- tibble(cell_type = coef_names) %>%
  mutate(cell_type = str_remove(cell_type, "assignments"),
         cell_type = str_remove(cell_type, ":final_dtDT"))

sig_counts <- all_cell_types %>%
  left_join(sig_counts, by = "cell_type") %>%
  mutate(n_sig = replace_na(n_sig, 0))

p_ablation <- ggplot(sig_counts,
                     aes(x = reorder(cell_type, n_sig),
                         y = n_sig)) +
  geom_bar(stat = "identity", fill = "#c20000", width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(x = NULL,
       y = "Number of significant DE genes",
       title = "DT ablation-induced differential expression") +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(paste0(perturbs_folder, "DT_ablation_DE_gene_counts.png"),
       p_ablation, width = 4, height = 3, dpi = 300)

ggsave(paste0(perturbs_folder, "DT_ablation_DE_gene_counts.svg"),
       p_ablation, width = 4, height = 3)


stem_ablation <- all_results_ablation %>%
  filter(cell_type == "assignmentsStem_cells:final_dtDT",
         adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice(1:10)

top_10_genes <- stem_ablation$genes



logcpm <- cpm(dge, log = TRUE, prior.count = 1)

logcpm_top <- logcpm[top_10_genes, ]

rownames(meta_data) <- colnames(logcpm_top)

meta_plot <- meta_data %>%
  filter(assignments == "Stem_cells")

logcpm_top <- logcpm_top[, rownames(meta_plot)]
meta_plot <- meta_plot[colnames(logcpm_top), ]

logcpm_top <- logcpm_top[, meta_plot$assignments == "Stem_cells"]
meta_plot <- meta_plot[meta_plot$assignments == "Stem_cells", ]

plot_df <- as.data.frame(t(logcpm_top))
plot_df$sample <- rownames(plot_df)

plot_df <- plot_df %>%
  pivot_longer(cols = all_of(top_10_genes),
               names_to = "gene",
               values_to = "expression") %>%
  left_join(
    meta_plot %>%
      mutate(sample = rownames(meta_plot)) %>%
      select(sample, study_age, final_dt),
    by = "sample"
  )


plot_df$condition <- factor(plot_df$final_dt, levels = c("C", "DT"))


plot_one_gene_ablation <- function(gene_name, df) {
  
  gene_df <- df %>% 
    filter(gene == gene_name)
  
  control_means <- gene_df %>%
    filter(condition == "C") %>%
    group_by(study_age) %>%
    summarise(control_mean = mean(expression), .groups = "drop")
  
  gene_df <- gene_df %>%
    left_join(control_means, by = "study_age") %>%
    mutate(log2FC_centered = expression - control_mean)
  
  ggplot(gene_df,
         aes(x = study_age,
             y = log2FC_centered,
             fill = condition)) +
    
    stat_summary(fun = mean,
                 geom = "bar",
                 position = position_dodge(0.7),
                 width = 0.65,
                 color = "black") +
    
    geom_jitter(aes(color = condition),
                position = position_jitterdodge(
                  jitter.width = 0.15,
                  dodge.width = 0.7),
                size = 2,
                alpha = 0.9) +
    
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    scale_fill_manual(values = c("C" = "grey75",
                                 "DT" = "#ed4747")) +
    
    scale_color_manual(values = c("C" = "black",
                                  "DT" = "#7a0000")) +
    
    theme_classic(base_size = 14) +
    labs(title = gene_name,
         y = "log2 Fold Change (vs study control)",
         x = NULL,
         fill = NULL,
         color = NULL) +
    
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
}



output_dir <- paste0(perturbs_folder, "Stem_cell_DT_gene_plots/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (g in top_10_genes) {
  
  p <- plot_one_gene_ablation(g, plot_df)
  safe_name <- gsub("[^A-Za-z0-9_]", "_", g)
  
  ggsave(paste0(output_dir, safe_name, "_SC_DT.png"),
         p, width = 4, height = 5, dpi = 300)
  
  ggsave(paste0(output_dir, safe_name, "_SC_DT.svg"),
         p, width = 4, height = 5)
}














library(DESeq2)
library(Seurat)
library(tidyverse)
#update R

library(reticulate)
# Specify the path to the .h5ad file
h5ad_file <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/pdatas_2026_02_14.h5ad"

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

#only keep entries where Author is in
author_list = c("Ruggiero-Ruff et al. (2024)",
                "Miles et al. (2023)",
                "Qian et al. (2024)",
                "Guo et al. (2025)")

expression_table <- expression_table[,meta_data$Author %in% author_list]
meta_data <- meta_data[meta_data$Author %in% author_list,]


hfd = c("SRX21241800","SRX21241801", #Ruggiero-Ruff et al. (2024)
        "SRX21986106","SRX21986108","SRX21986113", "SRX21986114","SRX21986115","SRX21986116", #Miles et al. (2023)
        "SRX24151434","SRX24151435",#Qian et al. (2024)
        "SRX31166358","SRX31166359" #Guo et al. (2025)
)


meta_data$hfd <- ifelse(meta_data$SRA_ID %in% hfd, "HFD", "C")

#for #
#for SRX24151438, SRX24151439 set as IREPKO
irepko = c("SRX24151438", "SRX24151439")
#remove these
expression_table <- expression_table[,!meta_data$SRA_ID %in% irepko]
meta_data <- meta_data[!meta_data$SRA_ID %in% irepko,]

print(unique(meta_data$assignments[meta_data$hfd  == "HFD"]))
#keep only these in meta_data and exp table
assignments_to_keep <- unique(meta_data$assignments[meta_data$hfd == "HFD"])

table(meta_data$Author,meta_data$Comp_sex)
#Author should be concat Author and Comp_sex
meta_data$Author <- paste(meta_data$Author, meta_data$Comp_sex)
#find what assignments there are where prop1_mut is M
assignments_that_occur_twice = table(meta_data$assignments)
assignments_that_occur_twice = names(assignments_that_occur_twice[assignments_that_occur_twice > 6])
assignments_to_keep <-assignments_that_occur_twice


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

expression_table <- sweep(expression_table, 2, meta_data$psbulk_n_cells, "/")

expression_table <- expression_table * median(meta_data$psbulk_n_cells)

exp_table = expression_table
meta_data = meta_data

exp_table <- exp_table[,meta_data$assignments %in% assignments_to_keep]
meta_data <- meta_data[meta_data$assignments %in% assignments_to_keep,]

dge = DGEList(counts = exp_table, group = meta_data$assignments)

# Create the design matrix
design <- model.matrix(~0 + assignments + assignments:hfd, data = meta_data)
#change to make.names
colnames(design) <- make.names(colnames(design))

#remove genes expressed in less than 10% of the samples or with less than 1000 total counts
keep_genes <- filterByExpr(dge, design, min.total.count = 50,min.prop= 0.2)
sum(keep_genes)

dge <- dge[keep_genes, ]
gene_names <- rownames(dge)

dge <- calcNormFactors(dge, method = "TMM")


library(variancePartition)
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ 0 +  assignments + assignments:hfd + (1 | Author)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, meta_data, BPPARAM = param)

# fit dream model with contrasts
fit <- dream(vobjDream, form, meta_data)#, L)

head(coef(fit))

fit <- variancePartition::eBayes(fit)

head(coef(fit))

top_genes <- variancePartition::topTable(fit, coef="assignmentsStem_cells:hfdHFD", number = Inf, sort.by = "P",adjust.method="BH", lfc=0)
top_genes$genes <- rownames(top_genes)
print(head(top_genes,20))



coef_names <- c("assignmentsCorticotrophs:hfdHFD", 
                "assignmentsEndothelial_cells:hfdHFD",
                "assignmentsGonadotrophs:hfdHFD",
                "assignmentsImmune_cells:hfdHFD",
                "assignmentsLactotrophs:hfdHFD",
                "assignmentsMelanotrophs:hfdHFD",
                "assignmentsSomatotrophs:hfdHFD",
                "assignmentsStem_cells:hfdHFD",
                "assignmentsThyrotrophs:hfdHFD")
                
#save perturb results to csv
perturbs_folder = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/"
for (coef in coef_names){
  print(coef)
  top_genes <- topTable(fit, coef=coef, number = Inf, sort.by = "P",adjust.method="BH", lfc=0.25)
  top_genes$genes <- rownames(top_genes)
  print(head(top_genes,20))
  write.csv(top_genes, paste(perturbs_folder,coef,"_HFD_dream.csv",sep=""))
}


library(tidyverse)

perturbs_folder <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/"

coef_names <- c("assignmentsCorticotrophs:hfdHFD", 
                "assignmentsEndothelial_cells:hfdHFD",
                "assignmentsGonadotrophs:hfdHFD",
                "assignmentsImmune_cells:hfdHFD",
                "assignmentsLactotrophs:hfdHFD",
                "assignmentsMelanotrophs:hfdHFD",
                "assignmentsSomatotrophs:hfdHFD",
                "assignmentsStem_cells:hfdHFD",
                "assignmentsThyrotrophs:hfdHFD")

# Load all files into one dataframe
all_results <- map_df(coef_names, function(coef){
  
  file <- paste0(perturbs_folder, coef, "_HFD_dream.csv")
  
  read.csv(file) %>%
    mutate(cell_type = coef)
})

sig_results <- all_results %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(
    cell_type = str_remove(cell_type, "assignments"),
    cell_type = str_remove(cell_type, ":hfdHFD"),
    direction = ifelse(logFC > 0, "Upregulated", "Downregulated")
  )

sig_counts <- sig_results %>%
  count(cell_type, direction)

all_cell_types <- tibble(cell_type = coef_names) %>%
  mutate(cell_type = str_remove(cell_type, "assignments"),
         cell_type = str_remove(cell_type, ":hfdHFD"))

sig_counts <- all_cell_types %>%
  left_join(sig_counts, by = "cell_type") %>%
  mutate(n = replace_na(n, 0),
         direction = replace_na(direction, "Upregulated"))

#levels should upregulated then down
sig_counts$direction <- factor(sig_counts$direction, levels = c("Upregulated", "Downregulated"))

library(ggplot2)
p <- ggplot(sig_counts, 
            aes(x = reorder(cell_type, n, sum), 
                y = n, 
                fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#63b3ed",
                               "Downregulated" = "#ffa500")) +
  theme_classic(base_size = 14) +
  labs(x = NULL,
       y = "Number of significant DE genes",
       fill = NULL,
       title = "HFD-induced DEG") +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p
#ggsave to perturbs
ggsave(filename = paste0(perturbs_folder, "HFD_DE_gene_counts.png"), plot = p, width = 5.5, height = 3, dpi = 300)
#svg
ggsave(filename = paste0(perturbs_folder, "HFD_DE_gene_counts.svg"), plot = p, width = 5.5, height = 3)




# Stem cell subset
stem_top <- all_results[all_results$cell_type == "assignmentsStem_cells:hfdHFD", ]

# Filter significant genes
top_10_genes <- stem_top[stem_top$adj.P.Val < 0.05, ]

#keep top 10
top_10_genes <- top_10_genes[1:10, ]

# Sort by adjusted p-value (most significant first)
top_10_genes <- top_10_genes[order(top_10_genes$adj.P.Val), ]

top_10_genes <- top_10_genes$genes


# log2 CPM (TMM normalized)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# Subset to top genes
logcpm_top <- logcpm[top_10_genes, ]

rownames(meta_data)<-colnames(logcpm_top)

meta_plot <- meta_data %>%
  filter(assignments == "Stem_cells")

# Make sure order matches dge columns
logcpm_top <- logcpm_top[,rownames(meta_plot )]
meta_plot <- meta_plot[colnames(logcpm_top), ]

# Keep only Stem cell columns
logcpm_top <- logcpm_top[, meta_plot$assignments == "Stem_cells"]
meta_plot <- meta_plot[meta_plot$assignments == "Stem_cells", ]

plot_df <- as.data.frame(t(logcpm_top))
plot_df$sample <- rownames(plot_df)

plot_df <- plot_df %>%
  pivot_longer(cols = all_of(top_10_genes),
               names_to = "gene",
               values_to = "expression") %>%
  left_join(
    meta_plot %>%
      mutate(sample = rownames(meta_plot)) %>%
      select(sample, Author, hfd),
    by = "sample"
  )

plot_df$hfd <- factor(plot_df$hfd, levels = c("C", "HFD"))

plot_df$Author <- factor(plot_df$Author,
                         levels = unique(plot_df$Author))

plot_df$gene <- factor(plot_df$gene,
                       levels = unique(plot_df$gene))

library(ggplot2)
library(ggplot2)
library(dplyr)
plot_one_gene <- function(gene_name, df) {
  
  gene_df <- df %>% 
    filter(gene == gene_name) %>%
    mutate(hfd = factor(hfd, levels = c("C", "HFD")))
  
  # Compute study-specific control means
  control_means <- gene_df %>%
    filter(hfd == "C") %>%
    group_by(Author) %>%
    summarise(control_mean = mean(expression), .groups = "drop")
  
  # Join and center expression
  gene_df <- gene_df %>%
    left_join(control_means, by = "Author") %>%
    mutate(log2FC_centered = expression - control_mean)
  
  p <- ggplot(gene_df,
              aes(x = Author,
                  y = log2FC_centered,
                  fill = hfd)) +
    
    # Mean bars
    stat_summary(fun = mean,
                 geom = "bar",
                 position = position_dodge(width = 0.7),
                 width = 0.65,
                 color = "black") +
    
    # Individual pseudobulk dots
    geom_jitter(aes(color = hfd),
                position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width = 0.7),
                size = 2,
                alpha = 0.9) +
    
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    
    scale_fill_manual(values = c("C" = "grey75",
                                 "HFD" = "#f79c54")) +
    
    scale_color_manual(values = c("C" = "black",
                                  "HFD" = "#b05b15")) +
    
    theme_classic(base_size = 14) +
    labs(title = gene_name,
         x = NULL,
         y = "log2 Fold Change (vs study control)",
         fill = NULL,
         color = NULL) +
    
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
  
  return(p)
}




gene_list <- unique(plot_df$gene)

gene_plots <- lapply(gene_list, function(g) {
  plot_one_gene(g, plot_df)
})



output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/Stem_cell_gene_plots/"

# Create folder if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_list <- unique(plot_df$gene)

for (g in gene_list) {
  
  p <- plot_one_gene(g, plot_df)
  
  # Clean filename (remove problematic characters)
  safe_name <- gsub("[^A-Za-z0-9_]", "_", g)
  
  # Save PNG (high-resolution)
  ggsave(
    filename = paste0(output_dir, safe_name, "_SC_HFD.png"),
    plot = p,
    width = 4,
    height = 5,
    dpi = 300
  )
  
  # Save SVG (vector format for Illustrator / Inkscape)
  ggsave(
    filename = paste0(output_dir, safe_name, "_SC_HFD.svg"),
    plot = p,
    width = 4,
    height = 5
  )
}






### Immune



stem_top <- all_results[all_results$cell_type == "assignmentsImmune_cells:hfdHFD", ]

# Filter significant genes
top_10_genes <- stem_top[stem_top$adj.P.Val < 0.05, ]

# Sort by adjusted p-value (most significant first)
top_10_genes <- top_10_genes[order(top_10_genes$adj.P.Val), ]

#keep top 10
top_10_genes <- top_10_genes[1:10, ]

top_10_genes <- top_10_genes$genes


# log2 CPM (TMM normalized)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# Subset to top genes
logcpm_top <- logcpm[top_10_genes, ]

rownames(meta_data)<-colnames(logcpm_top)

meta_plot <- meta_data %>%
  filter(assignments == "Immune_cells")

# Make sure order matches dge columns
logcpm_top <- logcpm_top[,rownames(meta_plot )]
meta_plot <- meta_plot[colnames(logcpm_top), ]

# Keep only Stem cell columns
logcpm_top <- logcpm_top[, meta_plot$assignments == "Immune_cells"]
meta_plot <- meta_plot[meta_plot$assignments == "Immune_cells", ]

plot_df <- as.data.frame(t(logcpm_top))
plot_df$sample <- rownames(plot_df)

plot_df <- plot_df %>%
  pivot_longer(cols = all_of(top_10_genes),
               names_to = "gene",
               values_to = "expression") %>%
  left_join(
    meta_plot %>%
      mutate(sample = rownames(meta_plot)) %>%
      select(sample, Author, hfd),
    by = "sample"
  )

plot_df$hfd <- factor(plot_df$hfd, levels = c("C", "HFD"))

plot_df$Author <- factor(plot_df$Author,
                         levels = unique(plot_df$Author))

plot_df$gene <- factor(plot_df$gene,
                       levels = unique(plot_df$gene))

library(ggplot2)
library(ggplot2)
library(dplyr)


gene_list <- unique(plot_df$gene)

gene_plots <- lapply(gene_list, function(g) {
  plot_one_gene(g, plot_df)
})


output_dir <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/perturbs/Stem_cell_gene_plots/"

# Create folder if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_list <- unique(plot_df$gene)

for (g in gene_list) {
  
  p <- plot_one_gene(g, plot_df)
  
  # Clean filename (remove problematic characters)
  safe_name <- gsub("[^A-Za-z0-9_]", "_", g)
  
  # Save PNG (high-resolution)
  ggsave(
    filename = paste0(output_dir, safe_name, "_Immune_HFD.png"),
    plot = p,
    width = 4,
    height = 5,
    dpi = 300
  )
  
  # Save SVG (vector format for Illustrator / Inkscape)
  ggsave(
    filename = paste0(output_dir, safe_name, "_Immune_HFD.svg"),
    plot = p,
    width = 4,
    height = 5
  )
}









