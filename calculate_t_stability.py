"""
calculate_t_stability.py - Advanced Cross-Domain Network Topology Calibration Pipeline.
Demonstrates the features of the pytstability engine across biology, linguistics, and culture.
"""

import os
import glob
import pandas as pd
import numpy as np

# Import the core calculator function from our package
from pytstability import calculate_densities

# Anchors all folders relative to this script's location
base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(base_dir, "data")
results_dir = os.path.join(base_dir, "results")
os.makedirs(results_dir, exist_ok=True)

print("\nStarting cross-domain network data pipeline")

# 1. Biological data example (Vertebrate Ecology)
print("\n" + "-" * 80)
print("\n Processing domain: vertebrate ecology (biology)\n")

bio_csvs = sorted(glob.glob(os.path.join(data_dir, "bio", "*.csv")))
bio_results = []

for csv_path in bio_csvs:
    label = str(os.path.basename(csv_path).split(".")[0])

    # Search for a matching tree file ending in .nex
    tree_path = csv_path.replace(".csv", "100trees.nex")
    if not os.path.exists(tree_path):
        tree_path = csv_path.replace(".csv", ".nex")

    local_csv_path = os.path.relpath(csv_path, base_dir)
    local_tree_path = os.path.relpath(tree_path, base_dir) if os.path.exists(tree_path) else "None"

    print(f"\n -> Processing file: {label}")
    print(f"    Spreadsheet: {local_csv_path}")
    print(f"    Tree file: {local_tree_path}")

    df_bio = pd.read_csv(csv_path)

    # get the metrics
    t_metrics = calculate_densities(
        data=df_bio,
        id_col='Species',
        tree=tree_path if os.path.exists(tree_path) else None,
        tree_type='fixed',
        struct_mode='tree',
        struct_depth='auto',
        verbose=True
    )

    bio_results.append({
        "Domain": "Biology",
        "Label": label,
        "Sample_Size_N": len(df_bio),
        "Spatial_Kappa": t_metrics["spatial_metrics"].get("kappa", "N/A"),
        "Spatial_P_Isolated": t_metrics["spatial_metrics"].get("P_isolated", "N/A"),
        "Spatial_S_Giant": t_metrics["spatial_metrics"].get("S_giant", "N/A"),
        "Structural_Kappa": t_metrics["structural_metrics"].get("kappa", "N/A"),
        "Structural_P_Isolated": t_metrics["structural_metrics"].get("P_isolated", "N/A"),
        "Structural_S_Giant": t_metrics["structural_metrics"].get("S_giant", "N/A"),
        "Topology_Stability_Index_TSI": t_metrics.get("topology_stability_index", "N/A"),
        "Critical_Horizon": t_metrics.get("critical_horizon", "N/A")
    })

if bio_results:
    bio_output_path = os.path.join(results_dir, "biology_connectivity_summary.csv")
    pd.DataFrame(bio_results).to_csv(bio_output_path, index=False)
    print(f"\nBiology results saved to: {os.path.relpath(bio_output_path, base_dir)}")


# 2. Linguistic data example (Linguistic Typology)
print("\n" + "-" * 80)
print("\n Processing domain: linguistic typology (linguistics)\n")

glottolog_path = os.path.join(data_dir, "tlu", "Glottolog_Languages.csv")
ling_results = []

if os.path.exists(glottolog_path):
    print(f" -> Loading language map files from: {os.path.relpath(glottolog_path, base_dir)}")
    gldf = pd.read_csv(glottolog_path)

    gldf = gldf.dropna(subset=['latitude', 'longitude'])
    gldf_shared = gldf[['glottocode', 'macroarea', 'Family_ID', 'longitude', 'latitude']].copy()
    gldf_shared['glottocode'] = gldf_shared['glottocode'].astype(str).str.strip().str.lower()

    feat_data_paths = sorted(glob.glob(os.path.join(data_dir, "tlu", "*", "*data.txt")))

    for feat_path in feat_data_paths:
        feat_id = str(os.path.basename(os.path.dirname(feat_path)))
        feat_dir = os.path.dirname(feat_path)
        tree_matches = glob.glob(os.path.join(feat_dir, "*trees.gz"))
        tree_path = tree_matches[0] if tree_matches else None

        print(f"\n -> Assessing language feature: {feat_id}")

        fdf = pd.read_csv(feat_path, delimiter="\t", header=None, names=["glottocode", "DV", "IV"])

        if len(fdf) < 10:
            print(f"    Skipping: not enough data rows (found {len(fdf)}, needs at least 10)")
            continue

        fdf['glottocode'] = fdf['glottocode'].astype(str).str.strip().str.lower()
        ling_df = pd.merge(gldf_shared, fdf, on='glottocode', how='inner')

        local_tree_print = os.path.relpath(tree_path, base_dir) if tree_path else "None"
        print(f"    Matched {len(ling_df)} languages on the map.")
        print(f"    Tree file found: {local_tree_print}")

        # get the metrics
        t_metrics = calculate_densities(
            data=ling_df,
            id_col='glottocode',
            tree=tree_path,
            tree_type='adaptive',
            struct_depth=8,
            coord_cols=['latitude', 'longitude'],
            spatial_km=500.0,
            spatial_cats='macroarea',
            struct_cats='Family_ID',
            verbose=True
        )

        ling_results.append({
            "Domain": "Linguistics",
            "Label": feat_id,
            "Sample_Size_N": len(ling_df),
            "Spatial_Kappa": t_metrics["spatial_metrics"].get("kappa", "N/A"),
            "Spatial_P_Isolated": t_metrics["spatial_metrics"].get("P_isolated", "N/A"),
            "Spatial_S_Giant": t_metrics["spatial_metrics"].get("S_giant", "N/A"),
            "Structural_Kappa": t_metrics["structural_metrics"].get("kappa", "N/A"),
            "Structural_P_Isolated": t_metrics["structural_metrics"].get("P_isolated", "N/A"),
            "Structural_S_Giant": t_metrics["structural_metrics"].get("S_giant", "N/A"),
            "Topology_Stability_Index_TSI": t_metrics.get("topology_stability_index", "N/A"),
            "Critical_Horizon": t_metrics.get("critical_horizon", "N/A")
        })

    if ling_results:
        ling_output_path = os.path.join(results_dir, "linguistics_connectivity_summary.csv")
        pd.DataFrame(ling_results).to_csv(ling_output_path, index=False)
        print(f"\nLinguistics results saved to: {os.path.relpath(ling_output_path, base_dir)}")
else:
    print(f"Cannot run linguistics loop: file missing at {os.path.relpath(glottolog_path, base_dir)}")


# 3. Cultural data example (Literary Themes)
print("\n" + "-" * 80)
print("\n Processing domain: cultural evolution (culture)\n")

etc_excel_paths = sorted(glob.glob(os.path.join(data_dir, "etc", "*.xlsx")))
etc_results = []

for excel_path in etc_excel_paths:
    label = str(os.path.basename(excel_path).split(".")[0])
    local_excel_path = os.path.relpath(excel_path, base_dir)
    print(f" -> Opening spreadsheet: {label}")
    print(f"    Path: {local_excel_path}")

    df_etc = pd.read_excel(excel_path)

    # get the metrics
    t_metrics = calculate_densities(
        data=df_etc,
        id_col='Name' if 'Name' in df_etc.columns else df_etc.columns[0],
        tree=None,
        tree_type='fixed',
        spatial_cats='Zone' if 'Zone' in df_etc.columns else None,
        struct_cats=['Century', 'Genre.1'] if 'Century' in df_etc.columns else None,
        verbose=True
    )

    etc_results.append({
        "Domain": "Culture",
        "Label": label,
        "Sample_Size_N": len(df_etc),
        "Spatial_Kappa": t_metrics["spatial_metrics"].get("kappa", "N/A"),
        "Spatial_P_Isolated": t_metrics["spatial_metrics"].get("P_isolated", "N/A"),
        "Spatial_S_Giant": t_metrics["spatial_metrics"].get("S_giant", "N/A"),
        "Structural_Kappa": t_metrics["structural_metrics"].get("kappa", "N/A"),
        "Structural_P_Isolated": t_metrics["structural_metrics"].get("P_isolated", "N/A"),
        "Structural_S_Giant": t_metrics["structural_metrics"].get("S_giant", "N/A"),
        "Topology_Stability_Index_TSI": t_metrics.get("topology_stability_index", "N/A"),
        "Critical_Horizon": t_metrics.get("critical_horizon", "N/A")
    })

if etc_results:
    etc_output_path = os.path.join(results_dir, "culture_connectivity_summary.csv")
    pd.DataFrame(etc_results).to_csv(etc_output_path, index=False)
    print(f"\nCulture results saved to: {os.path.relpath(etc_output_path, base_dir)}")

# =========================================================================
# 4. PIPELINE COMPLETE
# =========================================================================
print("\n" + "-" * 80)
print("Analysis complete. All multi-domain summaries have been successfully exported.\n")
