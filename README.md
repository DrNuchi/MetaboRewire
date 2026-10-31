# MetaboRewire <img src="MetaboRewire/man/figures/logo.png" align="right" height="139" alt="MetaboRewire logo" />

[![R-CMD-check](https://github.com/DrNuchi/MetaboRewire/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DrNuchi/MetaboRewire/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- After the first Zenodo release, add its DOI badge here -->

**Permutation-calibrated metabolite co-variation rewiring across multiple biological states (R).**

MetaboRewire is an R implementation and extension of
[PhenoRewire](https://github.com/ldellavedova/PhenoRewire) (Della Vedova *et al.*, bioRxiv 2026).
Its PhenoRewire core reproduces the Python implementation edge-for-edge (checked against
PhenoRewire v0.1.1 and v0.2.0; see `analysis/validation/`). It adds:

| Module | PhenoRewire | MetaboRewire |
|---|---|---|
| Feature selection | Spearman vs label, Monte-Carlo permutation, BH, bounded relaxation | same, plus **exact enumeration** of all labelings (deterministic) and optional limma selection |
| State networks | Spearman, adaptive abs(rho)/FDR, Louvain, centralities | identical |
| Rewiring | shared / state-only / sign-switch counts, node rewiring score; v0.2.0 adds a label-permutation null with fixed thresholds | exact label-permutation null that re-runs the whole adaptive procedure (global + node), leave-one-out edge support |
| Threshold-free | v0.2.0: connectivity-based selection at >= 14 samples per state | **differential connectivity (DC)** and **pathway cross-talk coherence** on the whole metabolome, usable at small n |
| More than two states | pairwise only | **three-state trajectory** (healthy -> disease -> treatment) with a regression-to-the-mean-corrected rescue null |
| Technical covariates | - | additive covariate (block/batch) adjustment |
| Several tissues | - | **paired dual-compartment coupling** and bipartite rewiring |

## Installation

Install the package directly from GitHub (the package lives in the `MetaboRewire/` subfolder):

```r
install.packages("remotes")
remotes::install_github("DrNuchi/MetaboRewire", subdir = "MetaboRewire")
# optional, for differential abundance and the study scripts:
install.packages(c("BiocManager", "vegan", "ggplot2", "patchwork", "ggrepel", "yaml"))
BiocManager::install("limma")
```

To reproduce the full study, download or clone the whole repository and follow the Quick start below.

## Data

The Metabolon data analysed in the study are not included in this repository.
Will be deposited in a repository soon but is available from the authors on request.

## Quick start (RStudio, macOS)

1. Install R (https://cloud.r-project.org) and RStudio.
2. Unzip and place the project next to your data, e.g.
   `/path/to/your/data` (code) and
   `/path/to/your/data` (Metabolon csv files).
3. Double-click **`MetaboRewire_project.Rproj`**, open **`RUN_ALL.R`**, check the SETTINGS block,
   click **Source**.

`RUN_ALL.R` installs missing packages, installs MetaboRewire from `MetaboRewire/`, runs the
unit tests and then the switched-on steps:

| Switch | Step | Script | Output |
|---|---|---|---|
| `run_main_analysis` | full study analysis (~15-25 min) | `analysis/01_run_eam_metformin.R` | `out_dir` |
| `run_simulation` | calibration / power simulation (~15 min on 4 cores) | `analysis/03_simulation_benchmark.R` | `out_dir/simulation` |
| `run_figures` | Figs 2-6, S1 (PDF + PNG) | `analysis/02_figures.R` | `figures/` |
| `run_supplementary` | Tables S1-S9 (CSV + LaTeX tabular) | `analysis/04_supplementary_tables.R` | `supplementary/` |
| `run_validation` | equivalence with Python PhenoRewire (needs Python, ~5 min) | `analysis/validation/*` | `analysis/validation/` |

Every run also writes `out_dir/provenance_<date>_<time>.txt` (configuration, steps, Python and
PhenoRewire versions, full R `sessionInfo()`), and `out_dir/sessionInfo.txt`.
`quick_test <- TRUE` does a 2-minute trial with few permutations (not for publication).

**Validation setup (once):** create a Python environment with PhenoRewire, then set
`python <- "/Users/<you>/phenorewire-env/bin/python"` and `run_validation <- TRUE`.

```bash
/opt/homebrew/bin/python3 -m venv ~/phenorewire-env
~/phenorewire-env/bin/python -m pip install "git+https://github.com/ldellavedova/PhenoRewire.git"
```

## Running from Terminal instead

```bash
cd MetaboRewire_project
Rscript -e 'install.packages(c("igraph","matrixStats","jsonlite","yaml","vegan","ggplot2","patchwork","ggrepel","testthat","BiocManager"), repos="https://cloud.r-project.org"); BiocManager::install("limma", update=FALSE)'
R CMD INSTALL MetaboRewire
cd analysis
Rscript 01_run_eam_metformin.R ../MetaboRewire/inst/config/eam_metformin.yaml
Rscript 03_simulation_benchmark.R 100 200 <out_dir>/simulation 4
Rscript 02_figures.R <out_dir> <out_dir>/simulation ../figures
Rscript 04_supplementary_tables.R <out_dir> validation <out_dir>/simulation ../supplementary
```

## Folder layout

```
MetaboRewire_project.Rproj     open this in RStudio
RUN_ALL.R                      one-click pipeline (RStudio: Source)
MetaboRewire/                  the R package
  R/                           source code (see "Code map" below)
  man/                         function documentation (?function_name in R)
  tests/testthat/              unit tests
  inst/config/eam_metformin.yaml   analysis settings for this study
analysis/                      analysis, figure, simulation and table scripts
  validation/                  R-vs-Python PhenoRewire equivalence check
  figures_src/                 Fig. 1 workflow diagram (Mermaid source + rendered PDF/PNG)
manuscript_template/           LaTeX templates (main text + supplement)
LICENSE                        MIT licence (with PhenoRewire notice)
CITATION.cff                   citation metadata ("Cite this repository" on GitHub)
.github/workflows/             automatic R CMD check on every push
```

## Code map (`MetaboRewire/R/`)

| File | Contents |
|---|---|
| `io_preprocess.R` | read Metabolon exports; presence filter, normalisation, log transform |
| `stats_selection.R` | rank-based association, exact / Monte-Carlo permutation, adaptive BH selection, Kruskal-Wallis |
| `networks.R` | Spearman networks, adaptive thresholds, topology metrics (PhenoRewire port) |
| `rewiring.R` | edge classes, rewiring summaries, triage score (PhenoRewire port) |
| `inference.R` | within-state standardisation, whole-procedure permutation test, jackknife support |
| `global_dc.R` | covariate adjustment, limma selection, whole-metabolome DC and pathway cross-talk |
| `trajectory.R` | three-state edge classes, rescue permutation tests, threshold-free rescue |
| `crosstalk_coupling.R` | pathway cross-talk enrichment, inter-organ coupling, bipartite rewiring |
| `pipeline.R` | one-call pairwise run (`run_phenorewire`), GraphML export |

## Outputs in `out_dir`

* `phenorewire/<variant>/<tissue>/<ref>_vs_<case>/` - PhenoRewire-format outputs: selected
  features, per-state `edges.csv` (with jackknife support), `priority_features.csv`,
  `rewiring/priority_rewired_nodes.csv` (rewiring score, DC, permutation p/q),
  `rewiring_permutation_global.csv`, `rewiring_network_phenotype.graphml` (Cytoscape), `summary.json`.
  Variants: `primary` (block-adjusted, PhenoRewire selection), `limma_selection`, `no_block_adjustment`.
* `whole_metabolome_dc/<tissue>/` - node DC with exact permutation p/q; pathway cross-talk.
* `three_state/<tissue>/` - three-state edge classes, rescue tests, threshold-free rescue, pathway-pair reversal.
* `heart_plasma_coupling/` - same-metabolite inter-organ coupling, coupling change, bipartite rewiring.
* `differential_abundance_<tissue>.csv`, `abundance_rescue_*.csv`, `permanova_*.csv`, `qc_samples_<tissue>.csv`.
* Summary tables: `phenorewire_pairwise_summary.csv`, `whole_metabolome_dc_summary.csv`,
  `three_state_*_summary.csv`.
* `metaborewire_results.rds` (all results as one R object, used by the figure script),
  `sessionInfo.txt`, `provenance_*.txt`.

## Using MetaboRewire on another dataset

```r
library(MetaboRewire)
ds <- read_metabolon_dataset("path/to/folder", suffix = "", tissue = "plasma")   # or build X/meta yourself
pp <- preprocess_matrix(ds$X, ds$detected, 0.5, groups = ds$meta$state)
X  <- adjust_covariate(pp$X, ds$meta$block, ds$meta$state)                       # optional
cfg <- yaml::read_yaml(system.file("config", "eam_metformin.yaml", package = "MetaboRewire"))$phenorewire
r  <- run_phenorewire(X, ds$meta, "Healthy", "EAM", ds$anno, ds$detected, cfg, "out/Healthy_vs_EAM")
cls <- setNames(ds$anno[rownames(X), "super_pathway"], rownames(X))
wm <- global_dc_test(X, ds$meta$state, c("Healthy", "EAM"), cls)
ts <- three_state_dc(X, ds$meta$state, c("Healthy", "EAM", "EAM_Met"), cls)
```

Any features x samples matrix works; annotation is only used to label and group nodes.

## Citation

If you use MetaboRewire, please cite [your paper / the Zenodo DOI] and PhenoRewire
(Della Vedova *et al.*, bioRxiv 2026, doi:10.64898/2026.08.29.748013).
GitHub's "Cite this repository" button (from `CITATION.cff`) gives ready-made formats.

## License

MIT. Algorithms in the PhenoRewire core are re-implemented from PhenoRewire (MIT,
(c) 2026 PhenoRewire contributors); please cite both.
