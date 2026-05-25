import os
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
import traceback
import sys
import os
# Add the epitome directory to the Python path
sys.path.append('/Users/k23030440/epitome_code/epitome/code')

# Import the required functions from the modules
from modules.heatmap import (
    process_heatmap_data,
    analyze_tf_cobinding,
    plot_heatmap
)

from modules.data_loader import load_heatmap_data

# Define the groupings with their descriptions
GROUPINGS = {
    "grouping_1_up": {"name": "Stem Cells vs Hormone-Producing Cells (UP)", "group_a": ["Stem_cells"], "group_b": ["Melanotrophs", "Corticotrophs", "Somatotrophs", "Lactotrophs", "Thyrotrophs", "Gonadotrophs"]},
    "grouping_1_down": {"name": "Stem Cells vs Hormone-Producing Cells (DOWN)", "group_a": ["Stem_cells"], "group_b": ["Melanotrophs", "Corticotrophs", "Somatotrophs", "Lactotrophs", "Thyrotrophs", "Gonadotrophs"]},
    "grouping_2_up": {"name": "Gonadotrophs vs Others (UP)", "group_a": ["Gonadotrophs"], "group_b": ["Stem_cells", "Melanotrophs", "Corticotrophs", "Somatotrophs", "Lactotrophs", "Thyrotrophs"]},
    "grouping_2_down": {"name": "Gonadotrophs vs Others (DOWN)", "group_a": ["Gonadotrophs"], "group_b": ["Stem_cells", "Melanotrophs", "Corticotrophs", "Somatotrophs", "Lactotrophs", "Thyrotrophs"]},
    "grouping_3_up": {"name": "TPIT-lineage vs Others (UP)", "group_a": ["Melanotrophs", "Corticotrophs"], "group_b": ["Stem_cells", "Gonadotrophs", "Somatotrophs", "Lactotrophs", "Thyrotrophs"]},
    "grouping_3_down": {"name": "TPIT-lineage vs Others (DOWN)", "group_a": ["Melanotrophs", "Corticotrophs"], "group_b": ["Stem_cells", "Gonadotrophs", "Somatotrophs", "Lactotrophs", "Thyrotrophs"]},
    "grouping_4_up": {"name": "Melanotrophs vs Corticotrophs (UP)", "group_a": ["Melanotrophs"], "group_b": ["Corticotrophs"]},
    "grouping_4_down": {"name": "Melanotrophs vs Corticotrophs (DOWN)", "group_a": ["Melanotrophs"], "group_b": ["Corticotrophs"]},
    "grouping_5_up": {"name": "PIT1-lineage vs Others (UP)", "group_a": ["Somatotrophs", "Lactotrophs", "Thyrotrophs"], "group_b": ["Stem_cells", "Gonadotrophs", "Melanotrophs", "Corticotrophs"]},
    "grouping_5_down": {"name": "PIT1-lineage vs Others (DOWN)", "group_a": ["Somatotrophs", "Lactotrophs", "Thyrotrophs"], "group_b": ["Stem_cells", "Gonadotrophs", "Melanotrophs", "Corticotrophs"]},
    "grouping_6_up": {"name": "Lactotrophs vs Other PIT1-lineage (UP)", "group_a": ["Lactotrophs"], "group_b": ["Somatotrophs", "Thyrotrophs"]},
    "grouping_6_down": {"name": "Lactotrophs vs Other PIT1-lineage (DOWN)", "group_a": ["Lactotrophs"], "group_b": ["Somatotrophs", "Thyrotrophs"]},
    "grouping_7_up": {"name": "Somatotrophs vs Other PIT1-lineage (UP)", "group_a": ["Somatotrophs"], "group_b": ["Lactotrophs", "Thyrotrophs"]},
    "grouping_7_down": {"name": "Somatotrophs vs Other PIT1-lineage (DOWN)", "group_a": ["Somatotrophs"], "group_b": ["Lactotrophs", "Thyrotrophs"]},
    "grouping_8_up": {"name": "Thyrotrophs vs Other PIT1-lineage (UP)", "group_a": ["Thyrotrophs"], "group_b": ["Lactotrophs", "Somatotrophs"]},
    "grouping_8_down": {"name": "Thyrotrophs vs Other PIT1-lineage (DOWN)", "group_a": ["Thyrotrophs"], "group_b": ["Lactotrophs", "Somatotrophs"]}
}

def extract_tf_hits(base_path, version="v_0.02", output_dir="results", generate_heatmaps=True):
    os.makedirs(output_dir, exist_ok=True)
    print(f"Loading data from {base_path}/heatmap/{version}...")
    
    try:
        motif_summary, coefs, rna_res, atac_res, mat, features, columns = load_heatmap_data(version=version)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # Use lists to accumulate to avoid overwriting or index issues
    master_all_list = []
    master_multimodal_list = []
    master_rna_list = []
    master_atac_list = []

    for grouping_key in GROUPINGS.keys():
        print(f"Processing {grouping_key}...")
        
        try:
            # Call process function
            # current_all_hits contains TFs even if motif_exists is False
            # plot_results contains TFs ONLY if motif_exists is True
            _, motifs, plot_results, current_all_hits = process_heatmap_data(
                motif_summary, coefs, rna_res, atac_res, mat, features, columns,
                grouping=grouping_key,
                AveExpr_threshold=1.5, mean_log2fc_threshold=0, fold_enrichment_threshold=0
            )

            # Define grouping/direction labels for this iteration
            direction = 'up' if grouping_key.endswith('_up') else 'down'
            grouping_dir = os.path.join(output_dir, grouping_key)
            os.makedirs(grouping_dir, exist_ok=True)

            # --- 1. Master List Accumulation ---
            # We use current_all_hits because plot_results filters out motif-less TFs
            iter_full = current_all_hits.copy()
            iter_full['grouping'] = grouping_key
            iter_full['direction'] = direction
            master_all_list.append(iter_full)
            iter_full.to_csv(os.path.join(grouping_dir, "all_hits.csv"), index=False)

            # --- 2. Multimodal Hits ---
            # Multimodal usually requires a motif, so we use plot_results
            m_subset = plot_results[plot_results['hit_type'] == 'multimodal'].copy()
            if not m_subset.empty:
                m_subset['grouping'] = grouping_key
                m_subset['direction'] = direction
                # Rename for standard output
                m_out = m_subset.rename(columns={
                    'p.adjust': 'atac_pvalue',
                    'geom_mean_adj_pval': 'rna_pvalue',
                    'fold.enrichment': 'atac_fold_enrichment',
                    'mean_log2fc': 'rna_log2fc'
                })
                master_multimodal_list.append(m_out)
                m_out.to_csv(os.path.join(grouping_dir, "multimodal_hits.csv"), index=False)

            # --- 3. RNA-only Hits ---
            # Use iter_full to catch TFs without motifs
            r_subset = iter_full[iter_full['hit_type'] == 'rna'].copy()
            if not r_subset.empty:
                r_out = r_subset.rename(columns={'geom_mean_adj_pval': 'rna_pvalue'})
                master_rna_list.append(r_out)
                r_out.to_csv(os.path.join(grouping_dir, "rna_hits.csv"), index=False)

            # --- 4. ATAC-only Hits ---
            a_subset = plot_results[plot_results['hit_type'] == 'atac'].copy()
            if not a_subset.empty:
                a_subset['grouping'] = grouping_key
                a_subset['direction'] = direction
                a_out = a_subset.rename(columns={
                    'p.adjust': 'atac_pvalue',
                    'fold.enrichment': 'atac_fold_enrichment'
                })
                master_atac_list.append(a_out)
                a_out.to_csv(os.path.join(grouping_dir, "atac_hits.csv"), index=False)

            # --- 5. Heatmaps ---
            if generate_heatmaps and len(motifs) > 0:
                from modules.heatmap import analyze_tf_cobinding
                res_df, fc_mat = analyze_tf_cobinding(_, motifs, return_matrix=True)
                res_df.to_csv(os.path.join(grouping_dir, "cobinding_stats.csv"), index=False)
                
                h_fig = plot_heatmap(res_df, motifs, plot_results, sig_threshold=0.05, fold_change_matrix=fc_mat)
                h_fig.savefig(os.path.join(grouping_dir, f"{grouping_key}_heatmap.png"), dpi=300, bbox_inches='tight')
                plt.close()

        except Exception:
            print(f"❌ Error processing {grouping_key}")
            traceback.print_exc()

    # --- FINAL COMBINATION ---
    print("\nSaving combined results...")
    if master_all_list:
        pd.concat(master_all_list, ignore_index=True).to_csv(os.path.join(output_dir, "all_hits_combined.csv"), index=False)
    if master_multimodal_list:
        pd.concat(master_multimodal_list, ignore_index=True).to_csv(os.path.join(output_dir, "multimodal_hits_combined.csv"), index=False)
    if master_rna_list:
        pd.concat(master_rna_list, ignore_index=True).to_csv(os.path.join(output_dir, "rna_hits_combined.csv"), index=False)
    if master_atac_list:
        pd.concat(master_atac_list, ignore_index=True).to_csv(os.path.join(output_dir, "atac_hits_combined.csv"), index=False)

    print("Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_path', type=str, required=True)
    parser.add_argument('--version', type=str, default='v_0.02')
    parser.add_argument('--output_dir', type=str, default='results')
    parser.add_argument('--no_heatmaps', action='store_true')
    args = parser.parse_args()
    
    extract_tf_hits(args.base_path, args.version, args.output_dir, not args.no_heatmaps)