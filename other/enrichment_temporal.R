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


#ageing res with temporal res
ageing_res = read.csv("/Users/k23030440/Library/CloudStorage/OneDrive-King'sCollegeLondon/PhD/Year_two/Aim 1/temporal_patterns/temporal_patterns.csv")


#perform the enrichment

library(enrichR)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(gridExtra)
#load ageing background 

# Base output folder for aging enrichment
aging_folder <- file.path(base_folder, "aging_temporal_enrichment")
dir.create(aging_folder, showWarnings = FALSE, recursive = TRUE)

# Set up enrichR (reusing the same databases)
setEnrichrSite("Enrichr")
dbs <- c("Reactome_Pathways_2024", "WikiPathways_2024_Mouse","GO_Biological_Process_2023")
cat("Using databases for aging analysis:", paste(dbs, collapse=", "), "\n")

# Get unique cell types
cell_types <- unique(ageing_res$cell_type)
cat(sprintf("Found %d unique cell types\n", length(cell_types)))

change_patterns = unique(ageing_res$change_pattern)
# Iterate through cell types and directions (up/down)
for(cell in cell_types) {
  
    
    #if not Stem_cells then skip
    if(cell != "Stem_cells") {
      next
    }
    
  for (dir in change_patterns) {
      cat(sprintf("Cell type: %s, Change pattern: %s\n", cell, dir))
      # Create analysis name and prepare directory
      analysis_name <- paste0("aging_", gsub(" ", "_", cell), "_", dir)
      output_folder <- file.path(aging_folder, analysis_name)
      dir.create(output_folder, showWarnings = FALSE)
      
      # Filter aging genes by cell type and direction
      filtered_genes <- ageing_res %>% 
        dplyr::filter(cell_type == cell & change_pattern == dir)
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
