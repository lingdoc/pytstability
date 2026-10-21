# pytstability

`pytstability` is a Python package that helps you figure out how well-connected your data is. It is built for researchers working in fields like biology, linguistics, and cultural anthropology.

Specifically, it looks at your dataset from two angles:
1. **Spatial Density:** How close your data points are to each other on a map.
2. **Structural Density:** How closely related your data points are historically or evolutionarily (like sharing branches on a family tree).

### Why use it?
When you feed a dataset into complex hierarchical models, the models assume your data points are well-connected. If your data is actually broken up into isolated "islands" with huge information gaps, those statistical models can fail—leading to massive errors, numbers that completely collapse, or topological mirages. This package flags those weak spots before you start modeling.

---

## Installation

You can install this package directly from GitHub using `pip`:

```bash
pip install git+https://github.com/lingdoc/pytstability
```

---

## Core Function Call & Outputs

To use the package, import the `calculate_densities` function:

```python
from pytstability import calculate_densities
```

### Unpacking Your Results
The engine returns **a single, comprehensive dictionary** containing all topological metrics and baseline diagnostic parameters. You catch the output using a single variable:

```python
results = calculate_densities(...)
```

The returned dictionary contains the following keys, which you can access directly:

*   **`results["spatial_metrics"]`** (dict): A sub-dictionary tracking your map connections, including overall density (`kappa`), the percentage of isolated data points (`P_isolated`), and the size of your largest single clump of data (`S_giant`).
*   **`results["structural_metrics"]`** (dict): A sub-dictionary tracking your tree or historical connections using the exact same metrics (`kappa`, `P_isolated`, `S_giant`).
*   **`results["topology_stability_index"]`** (float): The final **Topology Stability Index (TSI)**. This is your primary safety score telling you if the dataset is cohesive enough for standard modeling.
*   **`results["critical_horizon"]`** (float): Your baseline fallback safety limit (`gamma`), automatically calculated based on how geographically fragmented your data entries are.
*   **`results["spatial_components"]`** (int): The total number of separate geographic clusters discovered on your map.
*   **`results["structural_components"]`** (int): The total number of separate historical lineage islands found in your tree structure.

---

## Function Arguments

### Required Arguments
*   **`data`** (*pandas DataFrame*): The spreadsheet table containing your research data.
*   **`id_col`** (*string*): The column name that uniquely identifies each row (like `'Species'`, `'glottocode'`, or `'Language_ID'`).

### Optional Arguments
*   **`tree`** (*string*): The file path to your evolutionary or historical family tree file (supports `.nex`, `.csv`, and `.trees.gz`).
*   **`tree_type`** (*string*): Tells the code how to count steps down your tree branches.
    *   `"fixed"` (Default): Treats steps as flat integer counts. Best for standard biological species trees.
    *   `"adaptive"`: Converts steps into relative percentages. Perfect for language families where branches are uneven lengths.
*   **`struct_depth`** (*integer or `"auto"`*): How far back two points must share a branch to be considered "connected". **Defaults to `"auto"`, which lets the engine test your data and choose the best depth for you automatically.**
*   **`spatial_mode` / `struct_mode`** (*string*): Strategy for building the internal network grid (`'auto'`, `'coords'`, `'categorical'`, `'mixed'`, or `'none'`).
*   **`spatial_cats` / `struct_cats`** (*string or list*): Column headers used to group data by text matches (like matching a `'Zone'` or `'Family_ID'`).
*   **`coord_cols`** (*list of strings*): Your map coordinates, passed as a list, e.g., `['latitude', 'longitude']`.
*   **`spatial_km`** (*float*): The maximum distance in kilometers to consider two points geographically close. Defaults to `500.0` km.
*   **`verbose`** (*boolean*): Set to `True` to print progress logs and safety alerts straight to your terminal. Defaults to `True`.

---

## Feature Highlight: Let the Engine Autofit Your Settings

Picking a tree depth manually may bias your connectivity scores. If you set `struct_depth='auto'`, the engine takes over and tries to optimize the calculation settings for you, by finding the depth that leads to the best outcome:

*   **For Biology Datasets:** The code tests branch levels from 2 down to 10, running internal sanity checks to see which depth reveals the cleanest, most natural groupings.
*   **For Linguistics & Culture Datasets:** The code tests matching percentages from 50% to 98%. It cross-checks your historical tree against your geographic map points to find the setting where history and geography align most logically.

---

## Quick Code Examples

See the `calculate_k_densities.py` script for working code examples.

### 1. Vertebrate Ecology (Biology Tree with Autofit)
```python
import pandas as pd
from pytstability import calculate_densities

df = pd.read_csv("data/bio/mammals.csv")

# Letting struct_depth='auto' figure out the tree depth for you
results = calculate_densities(
    data=df,
    id_col="Species",
    tree="data/bio/mammals100trees.nex",  
    tree_type="fixed",              
    struct_depth="auto",
    verbose=True
)

# Extract metrics using the returned dictionary's keys
st_metrics = results["structural_metrics"]
tsi = results["topology_stability_index"]

print(f"Auto-selected Structural Density: {st_metrics['kappa']:.4f}")
print(f"Safety Stability Score (TSI): {tsi:.6f}")
```

### 2. Linguistic Typology (Combining Maps and Compressed Trees)
```python
import pandas as pd
from pytstability import calculate_densities

df = pd.read_csv("data/tlu/Glottolog_Languages.csv")

# The code automatically handles raw coordinates
results = calculate_densities(
    data=df,
    id_col="glottocode",
    tree="data/tlu/0015aKA/pruned_tree.trees.gz",
    tree_type="adaptive",
    struct_depth="auto",
    coord_cols=["latitude", "longitude"],
    spatial_cats="macroarea",
    struct_cats="Family_ID",
    spatial_km=500.0,
    verbose=True
)

print(f"Number of Spatial Clusters Found: {results['spatial_components']}")
print(f"Calculated Danger Horizon (Gamma): {results['critical_horizon']:.6f}")
```

### 3. Cultural Evolution (Category Groupings Only)
```python
import pandas as pd
from pytstability import calculate_densities

df = pd.read_excel("data/etc/Data_love_bywork_2019.xlsx")

# No tree file? No problem. It creates a network using your category text columns.
results = calculate_densities(
    data=df,
    id_col="Name",
    tree=None,
    spatial_cats="Zone",
    struct_cats=["Century", "Genre.1"],
    verbose=True
)

print(f"Dataset Topology Stability Index: {results['topology_stability_index']:.6f}")
```

---

## Reading Your Safety Reports

When `verbose=True`, the package prints out a quick diagnostic report about your data structure:

### PIPELINE CRASH ALERT
This means your dataset has broken apart into isolated islands. Your historical data is too fragmented, or your categories don't talk to each other. If you try to run traditional multi-level regression models on this data right now, your model's parameters are highly likely to collapse or blow up. You'll need to update your data, restructure your data partitions, or use specialized continuous models.

### ROBUST POSTERIOR SURFACE
This means your dataset passes the safety check. There is a healthy, continuous backbone of shared connections across your rows, making it perfectly safe to proceed with your standard regressions or MCMC parameter estimations.
