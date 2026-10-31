# =============================================================================
# MetaboRewire - run the whole study from RStudio
#
# HOW TO USE
#   1. Double-click MetaboRewire_project.Rproj (opens RStudio in this folder).
#   2. Open this file (RUN_ALL.R) and check the SETTINGS block below.
#   3. Click "Source" (top-right of the editor) - or press Cmd + Shift + S.
#
# Everything is printed in the Console. Total time: ~20-40 min for a full run,
# ~2 min with quick_test = TRUE (few permutations; for checking the setup only).
# =============================================================================

# ----------------------------- SETTINGS -------------------------------------
data_dir   <- "/Users/ovii/Desktop/metabolome/MetaboRewire"                     # folder with the Metabolon csv files
out_dir    <- "/Users/ovii/Desktop/metabolome/MetaboRewire/MetaboRewire_results"  # where results are written
quick_test <- FALSE   # TRUE = fast trial run with few permutations (results NOT for publication)
run_main_analysis <- TRUE
run_simulation    <- TRUE   # simulation benchmark (~15 min on 4 cores); written to <out_dir>/simulation
run_figures       <- TRUE
run_supplementary <- TRUE    # supplementary tables (CSV + LaTeX) into MetaboRewire_project/supplementary
run_validation    <- TRUE   # TRUE = compare with the original Python PhenoRewire (needs Python 3.10+; ~5 min)
python <- "/Users/ovii/phenorewire-env/bin/python"   # path to Python
cores             <- 4       # CPU cores for the simulation
# ----------------------------------------------------------------------------

project_dir <- getwd()
if (!file.exists(file.path(project_dir, "MetaboRewire", "DESCRIPTION")))
  stop("Working directory is not the project folder. Open MetaboRewire_project.Rproj ",
       "or run setwd('/path/to/MetaboRewire_project') first.")
if (!file.exists(file.path(data_dir, "batchnormimputed.csv")))
  stop("Cannot find batchnormimputed.csv in data_dir = ", data_dir, " - please fix data_dir in SETTINGS.")

# 1. Install missing packages ---------------------------------------------------
needed <- c("igraph", "matrixStats", "jsonlite", "yaml", "vegan", "ggplot2", "patchwork",
            "ggrepel", "testthat", "BiocManager")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma", update = FALSE, ask = FALSE)

# 2. Install MetaboRewire from this folder (re-installed every run, so code edits are picked up)
install.packages(file.path(project_dir, "MetaboRewire"), repos = NULL, type = "source", quiet = TRUE)
library(MetaboRewire)
message("MetaboRewire ", utils::packageVersion("MetaboRewire"), " installed. Running unit tests...")
testthat::test_dir(file.path(project_dir, "MetaboRewire", "tests", "testthat"), reporter = "summary")

# 3. Configuration -------------------------------------------------------------
cfg <- yaml::read_yaml(file.path(project_dir, "MetaboRewire", "inst", "config", "eam_metformin.yaml"))
cfg$data_dir <- data_dir
cfg$out_dir  <- if (quick_test) file.path(out_dir, "QUICK_TEST_not_for_publication") else out_dir
if (quick_test) {
  cfg$phenorewire$n_perm_rewiring <- 50
  cfg$whole_metabolome$n_perm <- 50
  cfg$three_state$n_perm <- 50
  cfg$three_state$n_perm_selection <- 200
  cfg$coupling$n_perm <- 50
}
cfg_file <- file.path(tempdir(), "metaborewire_config.yaml")
yaml::write_yaml(cfg, cfg_file)

run_script <- function(script, args, subdir = "analysis") {
  message("\n=== Running ", script, " ===")
  MR_ARGS <<- args
  on.exit(rm(MR_ARGS, envir = globalenv()), add = TRUE)
  old <- setwd(file.path(project_dir, subdir)); on.exit(setwd(old), add = TRUE)
  source(script, local = new.env(parent = globalenv()), echo = FALSE)
}

sim_dir <- file.path(cfg$out_dir, "simulation")
sim_for_figs <- function() {
  if (!file.exists(file.path(sim_dir, "simulation_summary.csv")))
    stop("No simulation results in ", sim_dir, ". Set run_simulation <- TRUE and run again.")
  sim_dir
}

# 4. Run --------------------------------------------------------------------------
t0 <- Sys.time()
if (run_main_analysis) run_script("01_run_eam_metformin.R", c(cfg_file, cfg$data_dir, cfg$out_dir))
if (run_simulation) run_script("03_simulation_benchmark.R",
                               c(if (quick_test) "5" else "100", if (quick_test) "50" else "200",
                                 sim_dir, as.character(cores)))
if (run_figures) run_script("02_figures.R", c(cfg$out_dir, sim_for_figs(),
                                               file.path(project_dir, if (quick_test) "figures_quick_test" else "figures")))
if (run_validation) {
  vdir <- file.path(project_dir, "analysis", "validation")
  py <- Sys.which(python); if (!nzchar(py)) py <- python
  ok_py <- tryCatch(system2(py, "--version", stdout = TRUE, stderr = TRUE), error = function(e) NULL)
  if (is.null(ok_py)) stop("Python not found. Install Python 3 (https://www.python.org/downloads/) ",
                           "and/or set `python` in SETTINGS to its full path.")
  message("\n=== Validation against Python PhenoRewire (", ok_py[1], ") ===")
  has_pr <- system2(py, c("-c", shQuote("import phenorewire")), stdout = FALSE, stderr = FALSE) == 0
  if (!has_pr) {
    message("Installing PhenoRewire from GitHub (needs git) ...")
    st <- system2(py, c("-m", "pip", "install", "--user",
                        "git+https://github.com/ldellavedova/PhenoRewire.git"))
    if (st != 0) stop("Could not install PhenoRewire. In Terminal run:\n  ", py,
                      " -m pip install git+https://github.com/ldellavedova/PhenoRewire.git")
  }
  run_script("make_validation_input.R", c(data_dir, "validation_input"), subdir = "analysis/validation")
  for (d in list.dirs(file.path(vdir, "validation_input"), recursive = FALSE)) {
    message("Python PhenoRewire: ", basename(d), " (a 'hard stop' error for some heart comparisons is expected)")
    system2(py, c(shQuote(file.path(vdir, "run_python_reference.py")), shQuote(d),
                  shQuote(file.path(vdir, "python_out", basename(d)))),
            stdout = file.path(vdir, paste0("python_", basename(d), ".log")),
            stderr = file.path(vdir, paste0("python_", basename(d), ".log")))
  }
  run_script("compare_r_python.R", c("validation_input", "python_out"), subdir = "analysis/validation")
  run_script("selection_seed_stability.R", "validation_input", subdir = "analysis/validation")
  message("Validation tables: ", file.path(vdir, "validation_R_vs_python.csv"), " and selection_seed_stability.csv")
}
if (run_supplementary)
  run_script("04_supplementary_tables.R", c(cfg$out_dir, file.path(project_dir, "analysis", "validation"),
                                            sim_for_figs(),
                                            file.path(project_dir, if (quick_test) "supplementary_quick_test" else "supplementary")))

# 5. Provenance: R session, package versions, settings, Python/PhenoRewire version -------
prov_file <- file.path(cfg$out_dir, paste0("provenance_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))
py_info <- tryCatch({
  py <- Sys.which(python); if (!nzchar(py)) py <- python
  c(system2(py, "--version", stdout = TRUE, stderr = TRUE),
    paste("PhenoRewire", system2(py, c("-c", shQuote("import importlib.metadata as m; print(m.version('phenorewire'))")),
                                 stdout = TRUE, stderr = TRUE)))
}, error = function(e) "Python / PhenoRewire: not available", warning = function(w) "Python / PhenoRewire: not available")
writeLines(c(paste("MetaboRewire run:", format(Sys.time())),
             paste("MetaboRewire version:", as.character(utils::packageVersion("MetaboRewire"))),
             paste("Steps run: main =", run_main_analysis, "| simulation =", run_simulation, "| figures =", run_figures,
                   "| supplementary =", run_supplementary, "| validation =", run_validation, "| quick_test =", quick_test),
             paste("data_dir:", data_dir), "", "--- Configuration used ---", readLines(cfg_file), "",
             "--- Python (validation) ---", py_info, "", "--- R sessionInfo() ---", capture.output(sessionInfo())),
           prov_file)
message("Provenance written: ", prov_file)
message("\nFinished in ", format(round(Sys.time() - t0, 1)), ".\nResults: ", cfg$out_dir,
        "\nFigures: ", file.path(project_dir, if (quick_test) "figures_quick_test" else "figures"))
