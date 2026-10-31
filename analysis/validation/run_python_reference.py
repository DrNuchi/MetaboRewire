"""Run the original PhenoRewire (Python) modules on a shared input and dump
intermediate results, bypassing the CLI hard stops (so that empty networks at
n = 6 do not abort the comparison)."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from phenorewire.preprocessing import prepare_matrices
from phenorewire.phenotype_selection import phenotype_selection
from phenorewire.networks import correlation_network, compute_network_metrics, compute_rewiring, summarize_rewiring

inp, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(exist_ok=True, parents=True)
df = pd.read_csv(inp / "feature_table.csv"); meta = pd.read_csv(inp / "metadata.csv")
gd = {"Reference": {"group_label": "reference"}, "Case": {"group_label": "case"}}
X, meta2, anno = prepare_matrices(df, meta, feature_id_col="feature_id", mz_col="mz", rt_col="rt",
    name_col="name", intensity_regex="^Sample_", meta_sample_col="sample_id", outdir=out,
    min_nonzero_intensity=0, min_presence_global=0.35, norm_method="median", log_transform=True,
    group_definition=gd)
y = (meta2["group"] == "Case").astype(int).values
_, sel, summ = phenotype_selection(X, y, anno, out, n_permutations=1000, fdr_alpha=0.05,
    adaptive_fdr=True, min_selected_features=5, fdr_alpha_ceiling=0.20, seed=42, min_features_hard_stop=1)
fs = sel.index.astype(str).tolist()
res = {"selection": summ}
edges = {}
for lab, mask in [("Reference", y == 0), ("Case", y == 1)]:
    G, e, th = correlation_network(X.loc[fs, X.columns[mask]], abs_r_thr=0.7, fdr_alpha=0.05,
        adaptive_thresholds=True, min_edges=3, max_density=0.25, floor_abs_r_thr=0.5, ceil_fdr_alpha=0.2)
    e.to_csv(out / f"edges_{lab}.csv", index=False); edges[lab] = e
    nodes, comm, ns = compute_network_metrics(sel, e, seed=42)
    nodes.to_csv(out / f"node_metrics_{lab}.csv")
    res[lab] = {"thresholds": th, "summary": ns.iloc[0].to_dict()}
sh, oa, ob, union = compute_rewiring(edges["Reference"], edges["Case"], label_a="Reference", label_b="Case", sign_switch_min_r=0.5)
rs, rn = summarize_rewiring(union, label_a="Reference", label_b="Case")
rs.to_csv(out / "rewiring_summary.csv", index=False); rn.to_csv(out / "rewiring_nodes.csv", index=False)
X.to_csv(out / "X_processed.csv")
json.dump(res, open(out / "summary.json", "w"), indent=1, default=str)
print(json.dumps(res, default=str)[:1500])
