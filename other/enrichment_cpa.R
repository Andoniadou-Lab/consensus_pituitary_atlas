# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

#update enrichR
#BiocManager::install("enrichR", ask = FALSE, update = TRUE)
library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(gridExtra)
library(purrr)

#load paste0("/Users/k23030440/epitome_code/epitome/markers/v_0.01/grouping_lineage_markers.csv"
markers = read.csv(file = "/Users/k23030440/epitome_code/epitome/data/markers/v_0.01/grouping_lineage_markers.csv", header = TRUE)

markers


# Base output folder
base_folder <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/enrichment"
dir.create(base_folder, showWarnings = FALSE, recursive = TRUE)

# Set up enrichR
setEnrichrSite("Enrichr")  # Use the main Enrichr site
#print possible databases
enrichr_dbs <- listEnrichrDbs()
dbs <- c("Reactome_Pathways_2024", "WikiPathways_2024_Mouse","GO_Biological_Process_2023")
cat("Using databases:", paste(dbs, collapse=", "), "\n")

# Base output folder
base_folder <- "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/enrichment"
dir.create(base_folder, showWarnings = FALSE, recursive = TRUE)

# Process and filter results function 
process_results <- function(results, database_name, pval_cutoff = 0.05) {
  if(is.null(results[[database_name]]) || nrow(results[[database_name]]) == 0) {
    return(NULL)
  }
  
  df <- results[[database_name]] %>%
    dplyr::filter(Adjusted.P.value < pval_cutoff) %>%
    dplyr::mutate(
      Database = database_name,
      logP = -log10(Adjusted.P.value),
      Term = str_replace_all(Term, "\\s*\\(.*?\\)\\s*", ""),  # Remove parenthetical content
      Term = str_replace_all(Term, "Homo sapiens|Mus musculus", ""),  # Remove species names
      Term = str_trim(Term)
    ) %>%
    arrange(Adjusted.P.value)
  
  return(df)
}

#load coef
coefs_table <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/coef.csv")
coefs_table
background = coefs_table$X

# Iterate through all groupings and directions
for(group_num in 1:8) {
  for(dir in c("up", "down")) {
    # Create analysis name and prepare directory
    analysis_name <- paste0("grouping_", group_num, "_", dir)
    output_folder <- file.path(base_folder, analysis_name)
    dir.create(output_folder, showWarnings = FALSE)
    
    # Filter directly from markers (which should be in the global environment)
    group_var <- paste0("grouping_", group_num)
    filtered_markers <- markers %>% 
      dplyr::filter(grouping == group_var & direction == dir)
    genes <- filtered_markers$gene
    
    # Print number of genes for verification
    cat(sprintf("\n\n===== Analyzing %s =====\n", analysis_name))
    cat(sprintf("Number of genes: %d\n", length(genes)))
    
    # Skip if no genes found
    if(length(genes) == 0) {
      cat(sprintf("No genes found for %s. Skipping...\n", analysis_name))
      next
    }
    
    # Perform enrichment analysis
    enrichment_results <- enrichr(genes, dbs, background = background, include_overlap=TRUE)
    
    # Process results from both databases
    reactome_results <- process_results(enrichment_results, "Reactome_Pathways_2024")
    wiki_results <- process_results(enrichment_results, "WikiPathways_2024_Mouse")
    go_results <- process_results(enrichment_results, "GO_Biological_Process_2023")
    
    # Combine results
    combined_results <- bind_rows(reactome_results, wiki_results, go_results)
    
    #remove terms with more than 200 genes or less than 4 genes - use "Genes" column
    combined_results <- combined_results %>%
      dplyr::mutate(n_genes = lengths(strsplit(Genes, ";"))) %>%
      dplyr::filter(n_genes <= 200 & n_genes >= 4)
    
    # Save to CSV in the group-specific folder
    results_file <- file.path(output_folder, "gene_enrichment_results.csv")
    write.csv(combined_results, results_file, row.names = FALSE)
    cat(sprintf("Full results saved to: %s\n", results_file))
    
    # Create visualization if results exist
    if(!is.null(combined_results) && nrow(combined_results) > 0) {
      
      combined_results <- combined_results %>%
        dplyr::mutate(genes_in_base_set = as.numeric(str_extract(`Overlap`, "(?<=/)[0-9]+"))) %>%
        dplyr::filter(genes_in_base_set <= 200)
      
      
      # Keep top 10 from each Database
      vis_results <- combined_results %>%
        group_by(Database) %>%
        slice_min(Adjusted.P.value, n = 5) %>%
        ungroup()
    
      
      # Reorder factors for plotting
      vis_results <- vis_results %>%
        dplyr::mutate(Term = factor(Term, levels = rev(unique(Term))))

      
      # Create plot
      # Define blue color palette
      blue_colors <- c("#ADD8E6", "#00008B", "#4169E1")
      names(blue_colors) <- unique(vis_results$Database)
      
      # Create plot with term truncation for display only, while preserving unique y positions
      p <- ggplot(vis_results, aes(x = logP, y = Term, fill = Database)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = blue_colors) +
        scale_y_discrete(labels = function(x) {
          # Create shortened labels for display only
          sapply(x, function(term) {
            if(nchar(term) > 40) {
              paste0(substr(term, 1, 40), "...")
            } else {
              term
            }
          })
        }) +
        labs(
          title = paste0("Pathway Enrichment Analysis: ", analysis_name),
          subtitle = sprintf("Top pathways from each database\n (p < 0.05) from %d input genes", 
                             length(genes)),
          x = expression(-log[10](adjusted~p-value)),
          y = NULL
        ) +
        theme_minimal() +
        theme(
          text = element_text(family = "sans", size = 12, color = "black"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 12, face = "italic"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
        )
      
      # Save plot as high-resolution PDF and PNG
      plot_svg <- file.path(output_folder, "pathway_enrichment.svg")
      plot_png <- file.path(output_folder, "pathway_enrichment.png")
      
      ggsave(plot_svg, p, width = 10, height = 8, units = "in", dpi = 300)
      ggsave(plot_png, p, width = 10, height = 8, units = "in", dpi = 300)
      
      cat(sprintf("Plots saved to:\n%s\n%s\n", plot_svg, plot_png))
      
      # Display plot
      print(p)
    } else {
      cat("No significant enrichment results found with p < 0.05\n")
    }
    
    # Print summary information
    cat("\nAnalysis Summary:\n")
    cat("------------------\n")
    cat(sprintf("Input genes: %d\n", length(genes)))
    cat(sprintf("Databases used: %s\n", paste(dbs, collapse=", ")))
    
    if(!is.null(reactome_results)) {
      cat(sprintf("Significant Reactome pathways: %d\n", nrow(reactome_results)))
    } else {
      cat("Significant Reactome pathways: 0\n")
    }
    
    if(!is.null(wiki_results)) {
      cat(sprintf("Significant WikiPathways: %d\n", nrow(wiki_results)))
    } else {
      cat("Significant WikiPathways: 0\n")
    }
    
    cat(sprintf("Results saved to: %s\n", results_file))
  }
}

# Create summary file with results from all analyses
cat("\n\nCreating overall summary...\n")
summary_file <- file.path(base_folder, "all_analyses_summary.csv")

# Prepare summary data structure
summary_rows <- list()

for(group_num in 1:8) {
  for(dir in c("up", "down")) {
    analysis_name <- paste0("grouping_", group_num, "_", dir)
    results_file <- file.path(base_folder, analysis_name, "gene_enrichment_results.csv")
    
    # Try to read results file if it exists
    if(file.exists(results_file)) {
      results <- read.csv(results_file)
      
      # Get gene count
      group_var <- paste0("grouping_", group_num)
      filtered_markers <- markers %>% 
        dplyr::filter(grouping == group_var & direction == dir)
      gene_count <- length(filtered_markers$gene)
      
      # Get top pathway if available
      if(nrow(results) > 0) {
        top_result <- results %>% arrange(Adjusted.P.value) %>% slice(1)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          analysis = analysis_name,
          input_genes = gene_count,
          significant_pathways = nrow(results),
          top_pathway = top_result$Term[1],
          top_pvalue = top_result$Adjusted.P.value[1],
          database = top_result$Database[1]
        )
      } else {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          analysis = analysis_name,
          input_genes = gene_count,
          significant_pathways = 0,
          top_pathway = NA,
          top_pvalue = NA,
          database = NA
        )
      }
    } else {
      # File doesn't exist, likely because no genes were found
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        analysis = analysis_name,
        input_genes = 0,
        significant_pathways = 0,
        top_pathway = NA,
        top_pvalue = NA,
        database = NA
      )
    }
  }
}

# Combine all summary rows and write to file
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, summary_file, row.names = FALSE)
cat(sprintf("Overall summary saved to: %s\n", summary_file))



#plotting a more biased list of genes for stem cells
results_sc <- read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/enrichment/grouping_1_up/gene_enrichment_results.csv")

#keep these results for the plot

results_sc = results_sc[results_sc$Term %in% c("Extracellular Matrix Organization",
                                          "Cell Junction Organization",
                                          "Cell-Cell Communication",
                                          
                                          "Focal Adhesion WP85",
                                          "Focal Adhesion PI3K Akt mTOR Signaling Pathway WP2841",
                                          "Oxidative Stress And Redox Pathway WP4466",
                                          
                                          "Positive Regulation Of Cell Population Proliferation",
                                          "Extracellular Matrix Organization",
                                          "Regulation Of Cell Migration",
                                          
                                          "ESC Pluripotency Pathways WP339",
                                          "Wnt Signaling Pathway WP539",
                                          "Ephrin Signaling",
                                          "Signaling by Hippo",
                                          "Hippo Signaling",
                                          "Canonical Wnt Signaling Pathway",
                                          "EGFR1 Signaling Pathway WP572",
                                          "Epidermal Growth Factor Receptor Signaling Pathway",
                                          "Non-Canonical Wnt Signaling Pathway",
                                          "Wnt Signaling Pathway",
                                          "TGF Beta Signaling Pathway WP113",
                                          "ErbB Signaling Pathway WP1261",
                                          "Laminin Interactions",
                                          "NOTCH3 Intracellular Domain Regulates Transcription",
                                          "Signaling by NOTCH3",
                                          "Notch Signaling Pathway",
                                          "ESC Pluripotency Pathways WP339",
                                          "FGFR1b Ligand Binding and Activation",
                                          "Signaling by TGFB Family Members",
                                          "YAP1- and WWTR1-stimulated Gene Expression"),]

vis_results <- results_sc %>%
  group_by(Term) %>%
  #take top 5 lowest pval
  top_n(1, -log10(Adjusted.P.value)) %>%
  ungroup()


# Reorder factor of Terms for plotting by sorting through pvals
vis_results <- vis_results %>%
  dplyr::mutate(Term = forcats::fct_reorder(Term, Adjusted.P.value, .desc = TRUE))


# Create plot with term truncation directly in the plot
# Define blue color palette
blue_colors <- c("#ADD8E6", "#00008B", "#4169E1")
names(blue_colors) <- unique(vis_results$Database)


# Create plot with term truncation for display only, while preserving unique y positions
p <- ggplot(vis_results, aes(x = logP, y = Term, fill = Database)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = blue_colors) +
  scale_y_discrete(labels = function(x) {
    # Create shortened labels for display only
    sapply(x, function(term) {
      if(nchar(term) > 35) {
        paste0(substr(term, 1, 35), "...")
      } else {
        term
      }
    })
  }) +
  labs(
    title = paste0("Pathway Enrichment Analysis: "),
    subtitle = ("Top pathways from each database"),
    x = expression(-log[10](adjusted~p-value)),
    y = NULL
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "sans", size = 12, color = "black"),
    plot.title = element_text(size = 0, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 12, face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

p

output_folder = "/Users/k23030440/Library/CloudStorage/OneDrive-King\'sCollegeLondon/PhD/Year_two/Aim\ 1/enrichment/"
# Save plot as high-resolution PDF and PNG
plot_svg <- file.path(output_folder, "pathway_enrichment.svg")
plot_png <- file.path(output_folder, "pathway_enrichment.png")

ggsave(plot_svg, p, width = 5, height = 6, units = "in", dpi = 300)
ggsave(plot_png, p, width = 5, height = 6, units = "in", dpi = 300)



#load ageing genes "/Users/k23030440/epitome_code/epitome/aging/",version,".csv"
ageing_res = read.csv("/Users/k23030440/epitome_code/epitome/data/aging/v_0.01/aging_genes.csv")
#keep only vals where adj.P.Val < 0.05
ageing_res = ageing_res[ageing_res$adj.P.Val < 0.05,]
#add col "direction based on logFC 
ageing_res$direction = ifelse(ageing_res$logFC > 0, "up", "down")

#perform the enrichment

library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(gridExtra)
#load ageing background 

background = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/enrichment/aging_enrichment/sc_aging_background_genes.csv"
background = read.csv(background)
background = background$x
# Base output folder for aging enrichment
aging_folder <- file.path(base_folder, "aging_enrichment")
dir.create(aging_folder, showWarnings = FALSE, recursive = TRUE)

# Set up enrichR (reusing the same databases)
setEnrichrSite("Enrichr")
dbs <- c("Reactome_Pathways_2024", "WikiPathways_2024_Mouse","GO_Biological_Process_2023")
cat("Using databases for aging analysis:", paste(dbs, collapse=", "), "\n")

# Get unique cell types
cell_types <- unique(ageing_res$cell_type)
cat(sprintf("Found %d unique cell types\n", length(cell_types)))

# Iterate through cell types and directions (up/down)
for(cell in cell_types) {
  #if not Stem_cells then skip
  if(cell != "Stem_cells") {
    next
  }
  
  for(dir in c("up", "down")) {
    # Create analysis name and prepare directory
    analysis_name <- paste0("aging_", gsub(" ", "_", cell), "_", dir)
    output_folder <- file.path(aging_folder, analysis_name)
    dir.create(output_folder, showWarnings = FALSE)
    
    # Filter aging genes by cell type and direction
    filtered_genes <- ageing_res %>% 
      dplyr::filter(cell_type == cell & direction == dir)
    genes <- filtered_genes$gene
    
    # Print number of genes for verification
    cat(sprintf("\n\n===== Analyzing %s =====\n", analysis_name))
    cat(sprintf("Number of aging genes: %d\n", length(genes)))
    
    # Skip if no genes found
    if(length(genes) < 10) {
      cat(sprintf("No genes found for %s. Skipping...\n", analysis_name))
      next
    }
    
    # Perform enrichment analysis
    enrichment_results <- enrichr(genes, dbs)#, background = background, include_overlap=TRUE)
    
    # Process results from both databases (reusing the process_results function)
    reactome_results <- process_results(enrichment_results, "Reactome_Pathways_2024")
    wiki_results <- process_results(enrichment_results, "WikiPathways_2024_Mouse")
    go_results <- process_results(enrichment_results, "GO_Biological_Process_2023")
    
    
    # Combine results
    combined_results <- bind_rows(reactome_results, wiki_results, go_results)
    
    
    # Save to CSV in the analysis-specific folder
    results_file <- file.path(output_folder, "gene_enrichment_results.csv")
    write.csv(combined_results, results_file, row.names = FALSE)
    cat(sprintf("Full results saved to: %s\n", results_file))
    
    # Create visualization if results exist
    if(!is.null(combined_results) && nrow(combined_results) > 0) {
      # Keep top 10 from each Database
      
      combined_results <- combined_results %>%
        dplyr::mutate(genes_in_base_set = as.numeric(str_extract(`Overlap`, "(?<=/)[0-9]+"))) %>%
        dplyr::filter(genes_in_base_set <= 300)
      
      #keep only those where Adjusted.P.value < 0.01
      vis_results <- combined_results[combined_results$Adjusted.P.value < 0.01,]
      
      
      vis_results <- vis_results %>%
        group_by(Database) %>%
        #take top 5 lowest pval
        top_n(5, -log10(Adjusted.P.value)) %>%
        ungroup()
        
      
      # Reorder factor of Terms for plotting by sorting through pvals
      vis_results <- vis_results %>%
        dplyr::mutate(Term = forcats::fct_reorder(Term, Adjusted.P.value, .desc = TRUE))
      
      
      # Create plot with term truncation directly in the plot
      # Define blue color palette
      blue_colors <- c("#ADD8E6", "#00008B", "#4169E1")
      names(blue_colors) <- unique(vis_results$Database)
      
      # Create plot with term truncation for display only, while preserving unique y positions
      p <- ggplot(vis_results, aes(x = logP, y = Term, fill = Database)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = blue_colors) +
        scale_y_discrete(labels = function(x) {
          # Create shortened labels for display only
          sapply(x, function(term) {
            if(nchar(term) > 40) {
              paste0(substr(term, 1, 40), "...")
            } else {
              term
            }
          })
        }) +
        labs(
          title = paste0("Pathway Enrichment Analysis: ", analysis_name),
          subtitle = sprintf("Top pathways from each database\n (p < 0.05) from %d input genes", 
                             length(genes)),
          x = expression(-log[10](adjusted~p-value)),
          y = NULL
        ) +
        theme_minimal() +
        theme(
          text = element_text(family = "sans", size = 12, color = "black"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = 9, color = "black"),
          axis.text.x = element_text(size = 9, color = "black"),
          axis.title.x = element_text(size = 12, face = "italic"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
        )
      
      # Save plot as high-resolution PDF and PNG
      plot_svg <- file.path(output_folder, "pathway_enrichment.svg")
      plot_png <- file.path(output_folder, "pathway_enrichment.png")
      
      ggsave(plot_svg, p, width = 5, height = 4.5, units = "in", dpi = 300)
      ggsave(plot_png, p, width = 5, height = 4.5, units = "in", dpi = 300)
      
      cat(sprintf("Plots saved to:\n%s\n%s\n", plot_svg, plot_png))
      
      # Display plot
      print(p)
    } else {
      cat("No significant enrichment results found with p < 0.05\n")
    }
    
    # Print summary information
    cat("\nAnalysis Summary:\n")
    cat("------------------\n")
    cat(sprintf("Input genes: %d\n", length(genes)))
    cat(sprintf("Databases used: %s\n", paste(dbs, collapse=", ")))
    
    if(!is.null(reactome_results)) {
      cat(sprintf("Significant Reactome pathways: %d\n", nrow(reactome_results)))
    } else {
      cat("Significant Reactome pathways: 0\n")
    }
    
    if(!is.null(wiki_results)) {
      cat(sprintf("Significant WikiPathways: %d\n", nrow(wiki_results)))
    } else {
      cat("Significant WikiPathways: 0\n")
    }
    
    cat(sprintf("Results saved to: %s\n", results_file))
  }
}

# Create a summary file for aging analysis
cat("\n\nCreating aging analysis summary...\n")
summary_file <- file.path(aging_folder, "aging_analysis_summary.csv")

# Prepare summary data structure
summary_rows <- list()

# Get unique cell types
cell_types <- unique(ageing_res$cell_type)

for(cell in cell_types) {
  for(dir in c("up", "down")) {
    analysis_name <- paste0("aging_", gsub(" ", "_", cell), "_", dir)
    results_file <- file.path(aging_folder, analysis_name, "gene_enrichment_results.csv")
    
    # Try to read results file if it exists
    if(file.exists(results_file)) {
      results <- read.csv(results_file)
      
      # Get gene count
      filtered_genes <- ageing_res %>% 
        dplyr::filter(cell_type == cell & direction == dir)
      gene_count <- nrow(filtered_genes)
      
      # Get top pathway if available
      if(nrow(results) > 0) {
        top_result <- results %>% arrange(Adjusted.P.value) %>% slice(1)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          analysis = analysis_name,
          cell_type = cell,
          direction = dir,
          input_genes = gene_count,
          significant_pathways = nrow(results),
          top_pathway = top_result$Term[1],
          top_pvalue = top_result$Adjusted.P.value[1],
          database = top_result$Database[1]
        )
      } else {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          analysis = analysis_name,
          cell_type = cell,
          direction = dir,
          input_genes = gene_count,
          significant_pathways = 0,
          top_pathway = NA,
          top_pvalue = NA,
          database = NA
        )
      }
    } else {
      # File doesn't exist, likely because no genes were found
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        analysis = analysis_name,
        cell_type = cell,
        direction = dir,
        input_genes = 0,
        significant_pathways = 0,
        top_pathway = NA,
        top_pvalue = NA,
        database = NA
      )
    }
  }
}

# Combine all summary rows and write to file
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, summary_file, row.names = FALSE)
cat(sprintf("Aging analysis summary saved to: %s\n", summary_file))

# Optional: Create an overall summary for aging genes (aggregating all cell types)
cat("\n\nCreating overall aging summary...\n")
overall_aging_summary_file <- file.path(aging_folder, "overall_aging_summary.csv")



#2 Months ageing

#perform the enrichment



#load ageing genes "/Users/k23030440/epitome_code/epitome/aging/",version,".csv"
ageing_res = read.csv("/Users/k23030440/epitome_code/epitome/data/aging/v_0.02/aging_genes_2_month.csv")
#keep only vals where adj.P.Val < 0.05
ageing_res = ageing_res[ageing_res$adj.P.Val < 0.05,]
#add col "direction based on logFC 
ageing_res$direction = ifelse(ageing_res$logFC > 0, "up", "down")

library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(gridExtra)
#load ageing background 

background = "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/enrichment/aging_enrichment/sc_aging_background_genes_2_month.csv"
background = read.csv(background)
background = background$x
# Base output folder for aging enrichment
aging_folder <- file.path(base_folder, "aging_enrichment")
dir.create(aging_folder, showWarnings = FALSE, recursive = TRUE)

# Set up enrichR (reusing the same databases)
setEnrichrSite("Enrichr")
dbs <- c("Reactome_Pathways_2024", "WikiPathways_2024_Mouse","GO_Biological_Process_2023")
cat("Using databases for aging analysis:", paste(dbs, collapse=", "), "\n")

# Get unique cell types
cell_types <- unique(ageing_res$cell_type)
cat(sprintf("Found %d unique cell types\n", length(cell_types)))

# Iterate through cell types and directions (up/down)
for(cell in cell_types) {
  #if not Stem_cells then skip
  if(cell != "Stem_cells") {
    next
  }
  
  for(dir in c("up", "down")) {
    # Create analysis name and prepare directory
    analysis_name <- paste0("aging_", gsub(" ", "_", cell), "_", dir)
    output_folder <- file.path(aging_folder, analysis_name)
    dir.create(output_folder, showWarnings = FALSE)
    
    # Filter aging genes by cell type and direction
    filtered_genes <- ageing_res %>% 
      dplyr::filter(cell_type == cell & direction == dir)
    genes <- filtered_genes$gene
    
    # Print number of genes for verification
    cat(sprintf("\n\n===== Analyzing %s =====\n", analysis_name))
    cat(sprintf("Number of aging genes: %d\n", length(genes)))
    
    # Skip if no genes found
    if(length(genes) < 10) {
      cat(sprintf("No genes found for %s. Skipping...\n", analysis_name))
      next
    }
    
    # Perform enrichment analysis
    enrichment_results <- enrichr(genes, dbs)#, background = background, include_overlap=TRUE)
    
    # Process results from both databases (reusing the process_results function)
    reactome_results <- process_results(enrichment_results, "Reactome_Pathways_2024")
    wiki_results <- process_results(enrichment_results, "WikiPathways_2024_Mouse")
    go_results <- process_results(enrichment_results, "GO_Biological_Process_2023")
    
    
    # Combine results
    combined_results <- bind_rows(reactome_results, wiki_results, go_results)
    
    
    # Save to CSV in the analysis-specific folder
    results_file <- file.path(output_folder, "gene_enrichment_results.csv")
    write.csv(combined_results, results_file, row.names = FALSE)
    cat(sprintf("Full results saved to: %s\n", results_file))
    
    # Create visualization if results exist
    if(!is.null(combined_results) && nrow(combined_results) > 0) {
      # Keep top 10 from each Database
      
      combined_results <- combined_results %>%
        dplyr::mutate(genes_in_base_set = as.numeric(str_extract(`Overlap`, "(?<=/)[0-9]+"))) %>%
        dplyr::filter(genes_in_base_set <= 300)
      
      #keep only those where Adjusted.P.value < 0.01
      vis_results <- combined_results[combined_results$Adjusted.P.value < 0.01,]
      
      
      vis_results <- vis_results %>%
        group_by(Database) %>%
        #take top 5 lowest pval
        top_n(5, -log10(Adjusted.P.value)) %>%
        ungroup()
      
      
      # Reorder factor of Terms for plotting by sorting through pvals
      vis_results <- vis_results %>%
        dplyr::mutate(Term = forcats::fct_reorder(Term, Adjusted.P.value, .desc = TRUE))
      
      
      # Create plot with term truncation directly in the plot
      # Define blue color palette
      blue_colors <- c("#ADD8E6", "#00008B", "#4169E1")
      names(blue_colors) <- unique(vis_results$Database)
      
      # Create plot with term truncation for display only, while preserving unique y positions
      p <- ggplot(vis_results, aes(x = logP, y = Term, fill = Database)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = blue_colors) +
        scale_y_discrete(labels = function(x) {
          # Create shortened labels for display only
          sapply(x, function(term) {
            if(nchar(term) > 40) {
              paste0(substr(term, 1, 40), "...")
            } else {
              term
            }
          })
        }) +
        labs(
          title = paste0("Pathway Enrichment Analysis: ", analysis_name),
          subtitle = sprintf("Top pathways from each database\n (p < 0.05) from %d input genes", 
                             length(genes)),
          x = expression(-log[10](adjusted~p-value)),
          y = NULL
        ) +
        theme_minimal() +
        theme(
          text = element_text(family = "sans", size = 12, color = "black"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = 9, color = "black"),
          axis.text.x = element_text(size = 9, color = "black"),
          axis.title.x = element_text(size = 12, face = "italic"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
        )
      
      # Save plot as high-resolution PDF and PNG
      plot_svg <- file.path(output_folder, "pathway_enrichment_2_month.svg")
      plot_png <- file.path(output_folder, "pathway_enrichment_2_month.png")
      
      ggsave(plot_svg, p, width = 5, height = 4.5, units = "in", dpi = 300)
      ggsave(plot_png, p, width = 5, height = 4.5, units = "in", dpi = 300)
      
      cat(sprintf("Plots saved to:\n%s\n%s\n", plot_svg, plot_png))
      
      # Display plot
      print(p)
    } else {
      cat("No significant enrichment results found with p < 0.05\n")
    }
    
    # Print summary information
    cat("\nAnalysis Summary:\n")
    cat("------------------\n")
    cat(sprintf("Input genes: %d\n", length(genes)))
    cat(sprintf("Databases used: %s\n", paste(dbs, collapse=", ")))
    
    if(!is.null(reactome_results)) {
      cat(sprintf("Significant Reactome pathways: %d\n", nrow(reactome_results)))
    } else {
      cat("Significant Reactome pathways: 0\n")
    }
    
    if(!is.null(wiki_results)) {
      cat(sprintf("Significant WikiPathways: %d\n", nrow(wiki_results)))
    } else {
      cat("Significant WikiPathways: 0\n")
    }
    
    cat(sprintf("Results saved to: %s\n", results_file))
  }
}

# Create a summary file for aging analysis
cat("\n\nCreating aging analysis summary...\n")
summary_file <- file.path(aging_folder, "aging_analysis_summary_2_month.csv")

# Prepare summary data structure
summary_rows <- list()

# Get unique cell types
cell_types <- unique(ageing_res$cell_type)

for(cell in cell_types) {
  for(dir in c("up", "down")) {
    analysis_name <- paste0("aging_", gsub(" ", "_", cell), "_", dir)
    results_file <- file.path(aging_folder, analysis_name, "gene_enrichment_results_2_month.csv")
    
    # Try to read results file if it exists
    if(file.exists(results_file)) {
      results <- read.csv(results_file)
      
      # Get gene count
      filtered_genes <- ageing_res %>% 
        dplyr::filter(cell_type == cell & direction == dir)
      gene_count <- nrow(filtered_genes)
      
      # Get top pathway if available
      if(nrow(results) > 0) {
        top_result <- results %>% arrange(Adjusted.P.value) %>% slice(1)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          analysis = analysis_name,
          cell_type = cell,
          direction = dir,
          input_genes = gene_count,
          significant_pathways = nrow(results),
          top_pathway = top_result$Term[1],
          top_pvalue = top_result$Adjusted.P.value[1],
          database = top_result$Database[1]
        )
      } else {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          analysis = analysis_name,
          cell_type = cell,
          direction = dir,
          input_genes = gene_count,
          significant_pathways = 0,
          top_pathway = NA,
          top_pvalue = NA,
          database = NA
        )
      }
    } else {
      # File doesn't exist, likely because no genes were found
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        analysis = analysis_name,
        cell_type = cell,
        direction = dir,
        input_genes = 0,
        significant_pathways = 0,
        top_pathway = NA,
        top_pvalue = NA,
        database = NA
      )
    }
  }
}

# Combine all summary rows and write to file
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, summary_file, row.names = FALSE)
cat(sprintf("Aging analysis summary saved to: %s\n", summary_file))

# Optional: Create an overall summary for aging genes (aggregating all cell types)
cat("\n\nCreating overall aging summary...\n")
overall_aging_summary_file <- file.path(aging_folder, "overall_aging_summary_2_month.csv")


























#these analyses/plots were not included in the publication

#load "/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes.csv"
sex_res = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/DE/sexually_dimorphic_genes.csv")
sex_res



# Base output folder for sex enrichment
sex_folder <- file.path(base_folder, "sex_enrichment")
dir.create(sex_folder, showWarnings = FALSE, recursive = TRUE)

# Set up enrichR (reusing the same databases)
setEnrichrSite("Enrichr")
dbs <- c("Reactome_Pathways_2024", "WikiPathways_2024_Mouse","GO_Biological_Process_2023")
cat("Using databases for sex analysis:", paste(dbs, collapse=", "), "\n")

# Get unique cell types
cell_types <- unique(sex_res$cell_type)
cat(sprintf("Found %d unique cell types\n", length(cell_types)))

# Iterate through cell types and directions (up/down)
for(cell in cell_types) {
  
  
  for(dir in c("male", "female")) {
    # Create analysis name and prepare directory
    analysis_name <- paste0("sex_", gsub(" ", "_", cell), "_", dir)
    output_folder <- file.path(sex_folder, analysis_name)
    dir.create(output_folder, showWarnings = FALSE)
    
    # Filter sex genes by cell type and direction
    filtered_genes <- sex_res %>% 
      dplyr::filter(cell_type == cell & sex == dir)
    genes <- filtered_genes$gene
    
    # Print number of genes for verification
    cat(sprintf("\n\n===== Analyzing %s =====\n", analysis_name))
    cat(sprintf("Number of sex genes: %d\n", length(genes)))
    
    # Skip if no genes found
    if(length(genes) < 10) {
      cat(sprintf("No genes found for %s. Skipping...\n", analysis_name))
      next
    }
    
    # Perform enrichment analysis
    enrichment_results <- enrichr(genes, dbs)
    
    # Process results from both databases (reusing the process_results function)
    reactome_results <- process_results(enrichment_results, "Reactome_Pathways_2024")
    wiki_results <- process_results(enrichment_results, "WikiPathways_2024_Mouse")
    go_results <- process_results(enrichment_results, "GO_Biological_Process_2023")
    
    
    # Combine results
    combined_results <- bind_rows(reactome_results, wiki_results, go_results)
    
    
    # Save to CSV in the analysis-specific folder
    results_file <- file.path(output_folder, "gene_enrichment_results.csv")
    write.csv(combined_results, results_file, row.names = FALSE)
    cat(sprintf("Full results saved to: %s\n", results_file))
    
    # Create visualization if results exist
    if(!is.null(combined_results) && nrow(combined_results) > 0) {
      
      #if less than 2 terms than skip
      if(nrow(combined_results) < 2) {
        cat("Less than 2 terms found. Skipping visualization...\n")
        next
      }
      # Keep top 10 from each Database
      
      combined_results <- combined_results %>%
        mutate(genes_in_base_set = as.numeric(str_extract(`Overlap`, "(?<=/)[0-9]+"))) %>%
        dplyr::filter(genes_in_base_set <= 300)
      
      #keep only those where Adjusted.P.value < 0.01
      vis_results <- combined_results[combined_results$Adjusted.P.value < 0.05,]
      
      
      vis_results <- vis_results %>%
        group_by(Database) %>%
        #take top 5 lowest pval
        top_n(5, -log10(Adjusted.P.value)) %>%
        ungroup()
        
      
      # Reorder factor of Terms for plotting by sorting through pvals
      vis_results <- vis_results %>%
        mutate(Term = forcats::fct_reorder(Term, Adjusted.P.value, .desc = FALSE))
      
      #if empty skip
      if(nrow(vis_results) == 0) {
        cat("No significant terms found. Skipping visualization...\n")
        next
      }
      
      # Create plot with term truncation directly in the plot
      # Define blue color palette
      blue_colors <- c("#ADD8E6", "#00008B", "#4169E1")
      names(blue_colors) <- unique(vis_results$Database)
      
      # Create plot with term truncation for display only, while preserving unique y positions
      p <- ggplot(vis_results, aes(x = logP, y = Term, fill = Database)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = blue_colors) +
        scale_y_discrete(labels = function(x) {
          # Create shortened labels for display only
          sapply(x, function(term) {
            if(nchar(term) > 40) {
              paste0(substr(term, 1, 40), "...")
            } else {
              term
            }
          })
        }) +
        labs(
          title = paste0("Pathway Enrichment Analysis: ", analysis_name),
          subtitle = sprintf("Top pathways from each database\n (p < 0.05) from %d input genes", 
                             length(genes)),
          x = expression(-log[10](adjusted~p-value)),
          y = NULL
        ) +
        theme_minimal() +
        theme(
          text = element_text(family = "sans", size = 12, color = "black"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 12, face = "italic"),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
        )
      
      # Save plot as high-resolution PDF and PNG
      plot_svg <- file.path(output_folder, "/pathway_enrichment.svg")
      plot_png <- file.path(output_folder, "/pathway_enrichment.png")
      
      ggsave(plot_svg, p, width = 5, height = 4, units = "in", dpi = 300)
      ggsave(plot_png, p, width = 5, height = 4, units = "in", dpi = 300)
      
      cat(sprintf("Plots saved to:\n%s\n%s\n", plot_svg, plot_png))
      
      # Display plot
      print(p)
    } else {
      cat("No significant enrichment results found with p < 0.05\n")
    }
    
    # Print summary information
    cat("\nAnalysis Summary:\n")
    cat("------------------\n")
    cat(sprintf("Input genes: %d\n", length(genes)))
    cat(sprintf("Databases used: %s\n", paste(dbs, collapse=", ")))
    
    if(!is.null(reactome_results)) {
      cat(sprintf("Significant Reactome pathways: %d\n", nrow(reactome_results)))
    } else {
      cat("Significant Reactome pathways: 0\n")
    }
    
    if(!is.null(wiki_results)) {
      cat(sprintf("Significant WikiPathways: %d\n", nrow(wiki_results)))
    } else {
      cat("Significant WikiPathways: 0\n")
    }
    
    cat(sprintf("Results saved to: %s\n", results_file))
  }
}

# Create a summary file for sex analysis
cat("\n\nCreating sex analysis summary...\n")
summary_file <- file.path(sex_folder, "sex_analysis_summary.csv")



filtered_genes <- sex_res %>% 
  dplyr::filter(cell_type == "Stem_cells" & sex == dir)
genes <- filtered_genes$gene
