# Seed dependence of Monte-Carlo permutation feature selection (PhenoRewire
# default, B = 1000) versus exact enumeration (MetaboRewire default) at n = 6/6.
suppressMessages(library(MetaboRewire))
# Arguments come from the command line (Rscript) or, in RStudio, from MR_ARGS set by RUN_ALL.R
args <- if (exists("MR_ARGS")) MR_ARGS else commandArgs(trailingOnly = TRUE)
inp_root <- if (length(args)) args[1] else "validation_input"
out <- list()
for (b in list.dirs(inp_root, recursive = FALSE, full.names = FALSE)) {
  ft <- read.csv(file.path(inp_root, b, "feature_table.csv"), check.names = FALSE)
  md <- read.csv(file.path(inp_root, b, "metadata.csv"))
  X <- as.matrix(ft[, grep("^Sample_", names(ft))]); rownames(X) <- ft$feature_id
  y <- as.integer(md$group_label == "case")
  pp <- preprocess_matrix(X, NULL, 0.35, norm = "median", transform = "log2p1")
  sets <- lapply(1:20, function(s) {
    r <- select_features(pp$X, y, "phenotype", n_perm = 1000, seed = s, min_features_hard_stop = 0)
    r$feature_id[r$selected] })
  ex <- select_features(pp$X, y, "phenotype", n_perm = "exact", min_features_hard_stop = 0)
  ex_set <- ex$feature_id[ex$selected]
  n <- sapply(sets, length)
  jac <- combn(20, 2, function(ij) { a <- sets[[ij[1]]]; bb <- sets[[ij[2]]]
    if (!length(union(a, bb))) 1 else length(intersect(a, bb)) / length(union(a, bb)) })
  out[[b]] <- data.frame(comparison = b, mc_min_selected = min(n), mc_median_selected = median(n),
                         mc_max_selected = max(n), mc_runs_with_zero = sum(n == 0),
                         mc_mean_pairwise_jaccard = round(mean(jac), 3), exact_selected = length(ex_set),
                         exact_applied_fdr = attr(ex, "meta")$applied_fdr_alpha)
}
res <- do.call(rbind, out)
write.csv(res, "selection_seed_stability.csv", row.names = FALSE)
print(res)
