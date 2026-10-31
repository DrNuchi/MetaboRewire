# Compare the MetaboRewire R port (PhenoRewire-compatibility settings) with the
# original Python PhenoRewire modules on identical inputs.
suppressMessages(library(MetaboRewire))
# Arguments come from the command line (Rscript) or, in RStudio, from MR_ARGS set by RUN_ALL.R
args <- if (exists("MR_ARGS")) MR_ARGS else commandArgs(trailingOnly = TRUE)
inp_root <- if (length(args)) args[1] else "validation_input"
py_root <- if (length(args) > 1) args[2] else "python_out"
rows <- list()
for (cmp in list.dirs(py_root, recursive = FALSE)) {
  b <- basename(cmp)
  # PhenoRewire <= 0.1.x writes phenotype_features_all.csv; 0.2.0 renamed it all_features_tested.csv
  feat_file <- c(file.path(cmp, "phenotype_features_all.csv"), file.path(cmp, "all_features_tested.csv"))
  feat_file <- feat_file[file.exists(feat_file)][1]
  if (!file.exists(file.path(cmp, "summary.json")) || is.na(feat_file)) {
    message("Skipping ", b, ": no Python output (see python_", b, ".log)"); next }
  ft <- read.csv(file.path(inp_root, b, "feature_table.csv"), check.names = FALSE)
  md <- read.csv(file.path(inp_root, b, "metadata.csv"))
  X <- as.matrix(ft[, grep("^Sample_", names(ft))]); rownames(X) <- ft$feature_id
  colnames(X) <- sub("^Sample_", "", colnames(X))
  y <- as.integer(md$group_label == "case")
  pp <- preprocess_matrix(X, NULL, min_presence_global = 0.35, norm = "median", transform = "log2p1")
  Xpy <- as.matrix(read.csv(file.path(cmp, "X_processed.csv"), row.names = 1, check.names = FALSE))
  colnames(Xpy) <- sub("^Sample_", "", colnames(Xpy))
  same_feats <- identical(sort(rownames(pp$X)), sort(rownames(Xpy)))
  max_dX <- max(abs(pp$X[rownames(Xpy), colnames(Xpy)] - Xpy))
  py_all <- read.csv(feat_file)
  py_all$feature_id <- as.character(py_all$feature_id)
  Xr <- pp$X[py_all$feature_id, ]
  r_mc <- perm_test_binary(Xr, y, n_perm = 1000, seed = 42)
  r_ex <- perm_test_binary(Xr, y, n_perm = "exact")
  max_dr <- max(abs(r_mc$r - py_all$r_pheno))
  cor_p_mc <- cor(r_mc$p, py_all$p_empirical); cor_p_ex <- cor(r_ex$p, py_all$p_empirical)
  sel_py <- py_all$feature_id[py_all$significant %in% c(TRUE, "True")]
  sel_r <- select_features(Xr, y, "phenotype", n_perm = 1000, fdr_alpha = 0.05, adaptive = TRUE,
                           min_selected = 5, alpha_ceiling = 0.2, min_features_hard_stop = 0)
  sel_rf <- sel_r$feature_id[sel_r$selected]
  jac <- length(intersect(sel_py, sel_rf)) / max(1, length(union(sel_py, sel_rf)))
  # networks on the Python-selected feature set: should match exactly
  p <- network_params()
  ekey <- function(e) if (nrow(e)) sort(paste(pmin(e$node_a, e$node_b), pmax(e$node_a, e$node_b))) else character()
  edge_match <- c(); r_diff <- c(); deg_diff <- c(); btw_diff <- c(); Q <- c()
  for (lab in c("Reference", "Case")) {
    mask <- if (lab == "Reference") y == 0 else y == 1
    nr <- correlation_network(pp$X[sel_py, mask, drop = FALSE], p)
    epy <- read.csv(file.path(cmp, paste0("edges_", lab, ".csv")), colClasses = c(node_a = "character", node_b = "character"))
    edge_match <- c(edge_match, identical(ekey(nr$edges), ekey(epy)))
    if (nrow(epy)) {
      k1 <- paste(pmin(nr$edges$node_a, nr$edges$node_b), pmax(nr$edges$node_a, nr$edges$node_b))
      k2 <- paste(pmin(epy$node_a, epy$node_b), pmax(epy$node_a, epy$node_b))
      r_diff <- c(r_diff, max(abs(nr$edges$r[match(k2, k1)] - epy$r)))
    }
    mr <- network_metrics(sel_py, nr$edges)
    mp <- read.csv(file.path(cmp, paste0("node_metrics_", lab, ".csv")))
    mp$feature_id <- as.character(mp[[1]])
    mm <- merge(mr$nodes, mp, by = "feature_id", suffixes = c("_r", "_py"))
    deg_diff <- c(deg_diff, max(abs(mm$degree_r - mm$degree_py)))
    btw_diff <- c(btw_diff, max(abs(mm$betweenness_centrality_r - mm$betweenness_centrality_py)))
    js <- jsonlite::fromJSON(file.path(cmp, "summary.json"))
    Q <- c(Q, sprintf("%.3f/%.3f", mr$summary$modularity, js[[lab]]$summary$modularity))
  }
  rows[[b]] <- data.frame(comparison = b, same_feature_set = same_feats, max_abs_diff_X = max_dX,
                          max_abs_diff_rho = max_dr, cor_p_R_MC_vs_py = round(cor_p_mc, 4),
                          cor_p_R_exact_vs_py = round(cor_p_ex, 4), n_sel_python = length(sel_py),
                          n_sel_R = length(sel_rf), jaccard_selection = round(jac, 3),
                          edges_identical = all(edge_match), max_abs_diff_edge_rho = if (length(r_diff)) max(r_diff) else 0,
                          max_abs_diff_degree = max(deg_diff), max_abs_diff_betweenness = max(btw_diff),
                          modularity_R_vs_py = paste(Q, collapse = "; "))
}
res <- do.call(rbind, rows)
if (is.null(res)) stop("No Python results found in ", py_root, ". Check the python_*.log files in this folder.")
write.csv(res, "validation_R_vs_python.csv", row.names = FALSE)
print(res)
