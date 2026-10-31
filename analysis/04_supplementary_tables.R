#!/usr/bin/env Rscript
# =============================================================================
# MetaboRewire - supplementary tables (CSV + LaTeX tabular snippets)
# Usage: Rscript 04_supplementary_tables.R <results_dir> <validation_dir> <simulation_dir> <out_dir>
#
# Writes one CSV per table (for Excel / the journal) and one .tex file per table
# containing only the tabular (no caption - write captions yourself).
# =============================================================================
args <- if (exists("MR_ARGS")) MR_ARGS else commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else "results/eam_metformin"
val_dir <- if (length(args) >= 2) args[2] else "validation"
sim_dir <- if (length(args) >= 3) args[3] else "../results/simulation"
out     <- if (length(args) >= 4) args[4] else "../supplementary"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

rd <- function(...) { f <- file.path(...); if (file.exists(f)) read.csv(f, check.names = FALSE) else NULL }
fmt <- function(df) {
  for (cn in names(df)) if (is.numeric(df[[cn]])) {
    x <- df[[cn]]
    df[[cn]] <- ifelse(is.na(x), "",
                ifelse(abs(x) > 0 & abs(x) < 1e-3, formatC(x, format = "e", digits = 1),
                ifelse(x == round(x), formatC(x, format = "d", big.mark = ""), formatC(x, format = "f", digits = 3))))
  }
  df
}
tex_escape <- function(x) gsub("([_%&#$])", "\\\\\\1", x)
write_tex <- function(df, file) {
  d <- fmt(df)
  al <- paste0("l", paste(ifelse(vapply(df, is.numeric, TRUE)[-1], "r", "l"), collapse = ""))
  lines <- c(sprintf("\\begin{tabular}{%s}", al), "\\toprule",
             paste0(paste(tex_escape(names(d)), collapse = " & "), " \\\\"), "\\midrule",
             apply(d, 1, function(r) paste0(paste(tex_escape(r), collapse = " & "), " \\\\")),
             "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file)
}
save_table <- function(df, name) {
  if (is.null(df) || !nrow(df)) { message("  - ", name, ": source not found, skipped"); return(invisible()) }
  utils::write.csv(df, file.path(out, paste0(name, ".csv")), row.names = FALSE)
  write_tex(df, file.path(out, paste0(name, ".tex")))
  message("  + ", name, " (", nrow(df), " rows)")
}

message("Writing supplementary tables to ", normalizePath(out))
# S1 - R port vs Python PhenoRewire
s1 <- rd(val_dir, "validation_R_vs_python.csv")
if (!is.null(s1)) s1 <- s1[, c("comparison", "same_feature_set", "max_abs_diff_X", "max_abs_diff_rho",
                               "cor_p_R_exact_vs_py", "n_sel_python", "n_sel_R", "edges_identical",
                               "max_abs_diff_degree", "max_abs_diff_betweenness", "modularity_R_vs_py")]
save_table(s1, "TableS1_R_vs_Python_equivalence")
# S2 - seed dependence of Monte-Carlo selection
save_table(rd(val_dir, "selection_seed_stability.csv"), "TableS2_selection_seed_stability")
# S3 - pairwise PhenoRewire runs with permutation calibration
s3 <- rd(res_dir, "phenorewire_pairwise_summary.csv")
if (!is.null(s3)) s3 <- s3[, c("tissue", "variant", "reference", "case", "n_selected", "applied_fdr",
                               "edges_ref", "edges_case", "state_specific", "state_specific_null_mean",
                               "p_state_specific", "rewiring_proportion", "rewiring_proportion_null_mean",
                               "DC", "DC_null_mean", "p_DC", "coherence_ref", "coherence_case",
                               "p_delta_coherence")]
save_table(s3, "TableS3_pairwise_rewiring")
# S4 - whole-metabolome differential connectivity
save_table(rd(res_dir, "whole_metabolome_dc_summary.csv"), "TableS4_whole_metabolome_DC")
# S5 - pathway cross-talk (all comparisons, both tissues)
ct <- do.call(rbind, lapply(c("plasma", "heart"), function(t) {
  fs <- list.files(file.path(res_dir, "whole_metabolome_dc", t), "^pathway_crosstalk_", full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) { d <- read.csv(f)
    data.frame(tissue = t, comparison = sub("pathway_crosstalk_(.*)\\.csv", "\\1", basename(f)), d) })) }))
if (!is.null(ct)) ct <- ct[order(ct$tissue, ct$comparison, ct$p), ]
save_table(ct, "TableS5_pathway_crosstalk")
# S6 - three-state analysis
s6a <- rd(res_dir, "three_state_threshold_free_summary.csv"); s6b <- rd(res_dir, "three_state_hard_summary.csv")
save_table(s6a, "TableS6a_three_state_threshold_free")
save_table(s6b, "TableS6b_three_state_hard_networks")
pr <- do.call(rbind, lapply(c("plasma", "heart"), function(t) { d <- rd(res_dir, "three_state", t, "pathway_crosstalk_rescue.csv")
  if (is.null(d)) NULL else data.frame(tissue = t, d) }))
save_table(pr, "TableS6c_pathway_pair_rescue")
# S7 - heart-plasma coupling
save_table(rd(res_dir, "heart_plasma_coupling", "same_metabolite_coupling.csv"), "TableS7a_heart_plasma_coupling")
save_table(rd(res_dir, "heart_plasma_coupling", "coupling_change_permutation.csv"), "TableS7b_coupling_change")
# S8 - simulation
save_table(rd(sim_dir, "simulation_summary.csv"), "TableS8_simulation_summary")
# S9 - global multivariate statistics
pm <- do.call(rbind, lapply(c("plasma", "heart"), function(t) rd(res_dir, paste0("permanova_pairwise_", t, ".csv"))))
save_table(pm, "TableS9_permanova_pairwise")
message("Done.")
