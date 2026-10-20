"""
pytstability/engine.py - Advanced Multi-Domain Network Topology Calibration Engine
Reviewed version with corrected spatial indexing and safe lineage translation fallbacks.
"""
import os
import io
import gzip
import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial.distance import pdist, squareform
from Bio import Phylo

def _compute_topological_metrics(A):
    """
    Computes advanced network topology metrics from an adjacency matrix.
    Returns: (edge_density, p_isolated, s_giant, num_components)
    """
    n = len(A)
    if n < 2:
        return 0.0, 1.0, 0.0, 1

    adj_no_self = A.copy()
    np.fill_diagonal(adj_no_self, 0.0)

    # 1. Edge Density (kappa)
    total_possible_edges = n * (n - 1)
    total_active_edges = np.sum(adj_no_self > 0)
    edge_density = float(total_active_edges / total_possible_edges) if total_possible_edges > 0 else 0.0

    # 2. Relative Proportion of Unlinked Nodes (P_isolated)
    node_degrees = np.sum(adj_no_self > 0, axis=1)
    num_isolated_nodes = np.sum(node_degrees == 0)
    p_isolated = float(num_isolated_nodes / n)

    # 3. Size of the Giant Connected Component and Component Counts
    try:
        sparse_graph = sp.csr_matrix(adj_no_self > 0)
        num_components, labels = sp.csgraph.connected_components(
            csgraph=sparse_graph, directed=False, return_labels=True
        )
        if num_components > 0:
            component_sizes = pd.Series(labels).value_counts()
            giant_component_size = component_sizes.max()
            s_giant = float(giant_component_size / n) if giant_component_size > 1 else 0.0
        else:
            s_giant = 0.0
            num_components = 1
    except Exception:
        s_giant = np.nan
        num_components = 1

    return edge_density, p_isolated, s_giant, num_components

def _audit_tree_segments(tree_path):
    """Extracts raw Newick path vectors to calculate scale-invariant branch variances."""
    open_func = gzip.open if str(tree_path).endswith('.gz') else open
    tree_string = ""

    with open_func(tree_path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_clean = line.strip()
            if line_clean.lower().startswith("tree ") or (line_clean.startswith("(") and line_clean.endswith(";")):
                tree_string = line_clean.split("=", 1)[1].strip() if "=" in line_clean else line_clean
                break

    if not tree_string:
        return np.nan

    tree_string = re.sub(r'\[.*?\]', '', tree_string)
    branch_matches = re.findall(r':([\d\.]+)', tree_string)
    lengths = np.array([float(b) for b in branch_matches])

    if len(lengths) == 0:
        return np.nan

    mean_segment = np.mean(lengths)
    var_segment = np.var(lengths)

    disconnection_ratio = float(var_segment / (mean_segment ** 2)) if mean_segment > 0 else 1.0
    return disconnection_ratio

def _parse_unified_tree(tree_path, df, id_col, tree_type='fixed', depth_threshold=3):
    """Reads a tree file, extracts the Newick structure, maps identifiers, and builds an adjacency matrix."""
    n = len(df)
    df_clean = df.copy()
    df_clean['match_key'] = df_clean[id_col].astype(str).str.strip().str.lower()
    species_list = df_clean['match_key'].tolist()

    open_func = gzip.open if str(tree_path).endswith('.gz') else open
    translate_map = {}
    tree_string = ""
    in_translate = False

    with open_func(tree_path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_clean = line.strip()
            if not line_clean:
                continue
            if line_clean.lower().startswith("translate"):
                in_translate = True
                continue
            if in_translate and line_clean == ";":
                in_translate = False
                continue
            if in_translate:
                parts = re.split(r'\s+', line_clean.rstrip(',').rstrip(';'))
                if len(parts) >= 2:
                    idx_token = parts[0].strip()
                    raw_label = parts[1].strip().replace("'", "").replace('"', '')
                    clean_label = raw_label.split('_')[0] if '_' in raw_label else raw_label
                    translate_map[idx_token] = clean_label.lower()
                continue
            if line_clean.lower().startswith("tree ") or (line_clean.startswith("(") and line_clean.endswith(";")):
                tree_string = line_clean.split("=", 1)[1].strip() if "=" in line_clean else line_clean
                break

    if not tree_string:
        raise ValueError(f"Could not find a valid Newick tree string in: {tree_path}")

    tree_string = re.sub(r'\[.*?\]', '', tree_string)
    lineage_paths_cache = {}

    try:
        target_tree = Phylo.read(io.StringIO(tree_string), "newick")
        for tip in target_tree.get_terminals():
            if not tip.name:
                continue
            tip_name = tip.name.strip().lower()
            resolved_key = translate_map.get(tip_name, tip_name).replace(' ', '_')
            lineage_paths_cache[resolved_key] = [target_tree.root] + target_tree.get_path(tip)
    except Exception:
        node_counter = 10000
        stack = [0]
        for i, char in enumerate(tree_string):
            if char == '(':
                stack.append(node_counter)
                node_counter += 1
            elif char == ')':
                if len(stack) > 1:
                    stack.pop()
            elif char not in [',', ';', ':', ' ']:
                match = re.match(r'^([\d\w\.\-_]+)', tree_string[i:])
                if match:
                    # Enforce strict sanitization strings on custom dictionary maps
                    token = match.group(1).strip().lower().replace("'", "").replace('"', '')
                    resolved_key = translate_map.get(token, token).replace(' ', '_')
                    lineage_paths_cache[resolved_key] = list(stack)

    A = np.zeros((n, n))
    mode = str(tree_type).strip().lower()
    cached_paths = [lineage_paths_cache.get(sp, lineage_paths_cache.get(sp.replace(' ', '_'), [0])) for sp in species_list]

    for i in range(n):
        path_i = cached_paths[i]
        len_i = len(path_i)
        for j in range(i + 1, n):
            path_j = cached_paths[j]
            len_j = len(path_j)
            min_length = min(len_i, len_j)
            shared_depth = sum(1 for k in range(min_length) if path_i[k] == path_j[k])

            if mode == 'adaptive':
                ratio_threshold = min(max(depth_threshold * 0.10, 0.05), 0.95)
                if (shared_depth / len_i if len_i > 0 else 0) > ratio_threshold and (shared_depth / len_j if len_j > 0 else 0) > ratio_threshold:
                    A[i, j] = A[j, i] = 1.0
            else:
                if shared_depth > depth_threshold:
                    A[i, j] = A[j, i] = 1.0
    return A

def calculate_densities(data, id_col, tree=None, tree_type='fixed',
                        spatial_mode='auto', struct_mode='auto',
                        spatial_cats=None, struct_cats=None, coord_cols=None,
                        spatial_km=500.0, struct_depth='auto', verbose=True):
    df = data.copy()
    n = len(df)

    spatial_metrics = {"kappa": np.nan, "P_isolated": np.nan, "S_giant": np.nan}
    struct_metrics = {"kappa": np.nan, "P_isolated": np.nan, "S_giant": np.nan}
    s_components = 1
    st_components = 1

    if n < 2:
        return spatial_metrics, struct_metrics, np.nan

    s_cats = [spatial_cats] if isinstance(spatial_cats, str) else list(spatial_cats) if spatial_cats else []
    st_cats = [struct_cats] if isinstance(struct_cats, str) else list(struct_cats) if struct_cats else []

    s_mode = str(spatial_mode).strip().lower()
    if s_mode == 'auto':
        if coord_cols and len(coord_cols) >= 2:
            s_mode = 'mixed' if s_cats else 'coords'
        else:
            s_mode = 'categorical' if s_cats else 'none'

    st_mode = str(struct_mode).strip().lower()
    if st_mode == 'auto':
        has_tree_file = tree is not None and os.path.exists(tree)
        if has_tree_file:
            st_mode = 'mixed' if st_cats else 'tree'
        else:
            st_mode = 'categorical' if st_cats else 'none'

    # =========================================================================
    # 3. Compile Spatial Adjacency Matrix Block
    # =========================================================================
    A_space = np.ones((n, n))
    calculated_space = False

    if s_mode in ['coords', 'mixed'] and coord_cols and len(coord_cols) >= 2:
        # implicit 2D array parsing to fix row metric stacks
        lat_rad = np.radians(df[[coord_cols[0]]].to_numpy().flatten().astype(float))
        lon_rad = np.radians(df[[coord_cols[1]]].to_numpy().flatten().astype(float))
        R = 6371.0
        coords = np.column_stack((R * np.cos(lat_rad) * np.cos(lon_rad), R * np.cos(lat_rad) * np.sin(lon_rad), R * np.sin(lat_rad)))
        A_space = (squareform(pdist(coords, metric='euclidean')) <= spatial_km).astype(float)
        calculated_space = True

    if s_mode in ['categorical', 'mixed'] and s_cats:
        A_space_cat = np.ones((n, n))
        valid_s_cols = [c for c in s_cats if c in df.columns]
        if valid_s_cols:
            for col in valid_s_cols:
                v = df[col].to_numpy()
                A_space_cat *= (v[:, None] == v[None, :]).astype(float)
            A_space = A_space * A_space_cat if calculated_space else A_space_cat
            calculated_space = True

    if calculated_space:
        k_s, p_s, s_g_s, s_components = _compute_topological_metrics(A_space)
        spatial_metrics = {"kappa": k_s, "P_isolated": p_s, "S_giant": s_g_s}

    # =========================================================================
    # 4. Compile Structural Adjacency Matrix Block (Parametric Framework)
    # =========================================================================
    A_tree = np.zeros((n, n))
    calculated_tree = False
    dr_base = 1.0

    if st_mode in ['tree', 'mixed'] and tree and os.path.exists(tree):
        try:
            dr_base = _audit_tree_segments(tree)
            if pd.isna(dr_base):
                dr_base = 1.0

            if str(struct_depth).strip().lower() == 'auto':
                best_entropy = -1.0
                best_param = 3
                best_A_tree = np.zeros((n, n))

                is_linguistic_or_cultural = (str(tree_type).strip().lower() == 'adaptive') or (s_mode in ['coords', 'mixed'])

                # ─── ROUTE 1: BIOLOGY DOMAINS ────────────────────────────────
                if not is_linguistic_or_cultural:
                    if verbose:
                        print(" Enforcing Fixed Absolute Depth Optimization (Biology Paradigm)...")
                    search_space = list(range(2, 11))
                    for test_d in search_space:
                        A_test = _parse_unified_tree(tree, df, id_col, 'fixed', test_d)
                        k_t, p_t, s_g_t, comp_t = _compute_topological_metrics(A_test)

                        s_g_clamped = max(s_g_t, 1e-6)
                        mlf_t = np.exp(-s_g_clamped / (k_t + 1e-5))
                        adjusted_components = (comp_t * mlf_t) / s_g_clamped
                        penalty_t = 1.0 + np.log(adjusted_components + 1.0)

                        test_tsi = (k_t * (1.0 - p_t)) / (dr_base * penalty_t)
                        topological_entropy = test_tsi * np.log(comp_t + 1.0) * (1.0 - p_t)

                        if topological_entropy > best_entropy:
                            best_entropy = topological_entropy
                            best_param = test_d
                            best_A_tree = A_test

                    struct_depth = best_param
                    tree_type = 'fixed'
                    if verbose:
                        print(f"   => Verified Optimal Fixed Depth: Level {struct_depth} (Entropy: {best_entropy:.4f})")

                # ─── ROUTE 2: LINGUISTICS DOMAINS ─────────────────────────────
                else:
                    if verbose:
                        print(" Enforcing Relative Ratio Threshold Optimization (Linguistics Paradigm)...")
                    search_space = [0.5, 0.7, 0.8, 0.85, 0.90, 0.93, 0.95, 0.96, 0.97, 0.98]
                    for test_r in search_space:
                        test_d_equiv = test_r / 0.10
                        A_test = _parse_unified_tree(tree, df, id_col, 'adaptive', test_d_equiv)
                        k_t, p_t, s_g_t, comp_t = _compute_topological_metrics(A_test)

                        s_g_clamped = max(s_g_t, 1e-6)
                        mlf_t = np.exp(-s_g_clamped / (k_t + 1e-5))
                        adjusted_components = (comp_t * mlf_t) / s_g_clamped
                        penalty_t = 1.0 + np.log(adjusted_components + 1.0)

                        test_tsi = (k_t * (1.0 - p_t)) / (dr_base * penalty_t)
                        base_entropy = test_tsi * np.log(comp_t + 1.0) * (1.0 - p_t)

                        if s_components > 1:
                            cross_domain_mismatch = abs(s_components - comp_t) / float(s_components)
                            cross_domain_penalty = np.exp(-cross_domain_mismatch * 3.5)
                            adjusted_entropy = base_entropy * cross_domain_penalty
                        else:
                            adjusted_entropy = base_entropy

                        if adjusted_entropy > best_entropy:
                            best_entropy = adjusted_entropy
                            best_param = test_r
                            best_A_tree = A_test

                    struct_depth = best_param
                    tree_type = 'adaptive'
                    if verbose:
                        print(f"   => Verified Optimal Adaptive Ratio: {best_param:.2f} (Adjusted Entropy: {best_entropy:.4f})")

                A_tree = best_A_tree
                calculated_tree = True
            else:
                A_tree = _parse_unified_tree(tree, df, id_col, tree_type, int(struct_depth))
                calculated_tree = True
        except Exception as e:
            if verbose: print(f"Warning: Tree parsing failed: {e}")

    A_family = np.zeros((n, n))
    calculated_family = False
    if st_mode in ['categorical', 'mixed'] and st_cats:
        valid_st_cols = [c for c in st_cats if c in df.columns]
        if valid_st_cols:
            for col in valid_st_cols:
                v = df[col].to_numpy()
                A_family += (v[:, None] == v[None, :]).astype(float)
            A_family = (A_family > 0).astype(float)
            calculated_family = True

    if calculated_tree or calculated_family:
        if calculated_tree and calculated_family:
            A_struct_final = ((A_tree + A_family) > 0).astype(float)
        elif calculated_tree:
            A_struct_final = A_tree
        else:
            A_struct_final = A_family

        k_t, p_t, s_g_t, st_components = _compute_topological_metrics(A_struct_final)
        struct_metrics = {"kappa": k_t, "P_isolated": p_t, "S_giant": s_g_t}

    # =========================================================================
    # 5. Combined TSI Formulation Layer & Dynamic Horizon
    # =========================================================================
    k_active = struct_metrics["kappa"] if not pd.isna(struct_metrics["kappa"]) else 0.0
    p_active = struct_metrics["P_isolated"] if not pd.isna(struct_metrics["P_isolated"]) else 0.0
    s_giant_active = struct_metrics["S_giant"] if not pd.isna(struct_metrics["S_giant"]) else 1.0

    s_giant_clamped = max(s_giant_active, 1e-6)
    manifold_linkage_factor = np.exp(-s_giant_clamped / (k_active + 1e-5))

    adjusted_components = (st_components * manifold_linkage_factor) / s_giant_clamped
    partition_penalty = 1.0 + np.log(adjusted_components + 1.0)

    dr_weighted = dr_base * partition_penalty
    topology_stability_index = float((k_active * (1.0 - p_active)) / dr_weighted)

    if verbose:
        print("")
        print("-- Data Topology Metrics --")
        print(f"Spatial Mode       : {s_mode} (Components: {s_components})")
        print(f"Structural Mode    : {st_mode} (Components: {st_components} | Optimized Depth: {struct_depth})")
        print(f"Spatial Topology   : κ = {spatial_metrics['kappa']:.4f} | P_isolated = {spatial_metrics['P_isolated']:.4f} | S_giant = {spatial_metrics['S_giant']:.4f}")
        print(f"Structural Topology: κ = {struct_metrics['kappa']:.4f} | P_isolated = {struct_metrics['P_isolated']:.4f} | S_giant = {struct_metrics['S_giant']:.4f}")
        print(f"Unified Diagnostic : Topology Stability Index (TSI) = {topology_stability_index:.6f}")
        print("")

    if not pd.isna(struct_metrics["kappa"]) and verbose:
        critical_horizon = 0.0050 * np.log(float(s_components) + 1.0)

        if topology_stability_index <= critical_horizon:
            print(f"* PIPELINE CRASH ALERT * Data saturation surface collapsed (Horizon: {critical_horizon:.5f}).")
            print(f"   Diagnostic Profile: TSI ({topology_stability_index:.5f}) resides inside the Critical Danger Zone.")
            print("   Required Action   : Topology-aware routing architectures or specialized GPGLMM wrappers are mandatory.")
        else:
            print(f"* PIPELINE STATUS *  : Robust posterior surface (Horizon: {critical_horizon:.5f}).")
            print("   Diagnostic Profile: TSI metrics confirm sufficient continuous shared variance backbones.")
        print("")

    return spatial_metrics, struct_metrics, topology_stability_index
