# ---- MetaboRewire - Association statistics and feature selection -------------
# Port of PhenoRewire stats.py / selection.py with two additions:
#  * exact enumeration of all label assignments for two-group designs when the
#    number of distinct assignments is small (here choose(12, 6) = 924), giving
#    exact rather than Monte-Carlo permutation p-values;
#  * a three-group Kruskal-Wallis permutation selector for multi-state designs.

.TIE_TOL <- 1e-12

#' Row-wise average ranks, centred and L2-normalised (so that their cross-product is Spearman rho)
#' @param X Numeric matrix, features x samples.
rank_rows <- function(X) {
  X <- as.matrix(X)
  R <- matrixStats::rowRanks(X, ties.method = "average")
  storage.mode(R) <- "double"
  R <- R - rowMeans(R)
  nr <- sqrt(rowSums(R^2))
  R[nr > 0, ] <- R[nr > 0, , drop = FALSE] / nr[nr > 0]
  R[nr == 0, ] <- 0
  R
}

.norm_vec <- function(y) { r <- rank(y, ties.method = "average"); r <- r - mean(r)
  if (sum(r^2) > 0) r / sqrt(sum(r^2)) else r }

#' Spearman correlation of every row of X with target y
#' @param X Numeric matrix, features x samples.
#' @param y Binary group indicator (0 = reference, 1 = case), one per sample.
spearman_assoc <- function(X, y) as.vector(rank_rows(X) %*% .norm_vec(y))

#' All (or a random sample of) label permutations of a two-group design.
#' Returns a logical matrix: n_perm x n_samples, TRUE = assigned to group 1.
#' @param y Binary group indicator (0 = reference, 1 = case), one per sample.
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param max_exact Maximum number of labelings for which exact enumeration is used.
#' @param seed Random seed (for reproducibility).
label_permutations <- function(y, n_perm = "exact", max_exact = 20000, seed = 42) {
  y <- as.logical(y); n <- length(y); k <- sum(y)
  n_tot <- choose(n, k)
  if (identical(n_perm, "exact") && n_tot <= max_exact) {
    cmb <- utils::combn(n, k)
    P <- matrix(FALSE, ncol(cmb), n)
    P[cbind(rep(seq_len(ncol(cmb)), each = k), as.vector(cmb))] <- TRUE
    attr(P, "exact") <- TRUE
    return(P)
  }
  if (identical(n_perm, "exact")) n_perm <- 2000
  set.seed(seed)
  P <- t(replicate(as.integer(n_perm), sample(y)))
  attr(P, "exact") <- FALSE
  P
}

#' Spearman association with a binary phenotype + permutation p-values
#'
#' Exact mode enumerates every assignment (observed included), p = #(|r*|>=|r|)/N.
#' Monte-Carlo mode follows PhenoRewire: p = (#(|r*|>=|r|) + 1) / (B + 1).
#' @param X Numeric matrix, features x samples.
#' @param y Binary group indicator (0 = reference, 1 = case), one per sample.
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
#' @param max_exact Maximum number of labelings for which exact enumeration is used.
perm_test_binary <- function(X, y, n_perm = "exact", seed = 42, max_exact = 20000) {
  Xn <- rank_rows(X)
  r_obs <- as.vector(Xn %*% .norm_vec(as.numeric(y)))
  P <- label_permutations(y, n_perm, max_exact, seed)
  Yn <- apply(P, 1, function(z) .norm_vec(as.numeric(z)))   # samples x perms
  Rp <- abs(Xn %*% Yn)
  cnt <- rowSums(Rp >= abs(r_obs) - .TIE_TOL)
  exact <- isTRUE(attr(P, "exact"))
  p <- if (exact) cnt / ncol(Yn) else (cnt + 1) / (ncol(Yn) + 1)
  p[!is.finite(r_obs)] <- 1
  list(r = r_obs, p = p, n_perm = ncol(Yn), exact = exact)
}

#' Spearman association with a numeric target (PhenoRewire temporal mode)
#' @param X Numeric matrix, features x samples.
#' @param t Numeric target (e.g. time), one per sample.
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
perm_test_numeric <- function(X, t, n_perm = 1000, seed = 42) {
  Xn <- rank_rows(X); tn <- .norm_vec(t)
  r_obs <- as.vector(Xn %*% tn)
  set.seed(seed)
  Tn <- replicate(n_perm, sample(tn))
  cnt <- rowSums(abs(Xn %*% Tn) >= abs(r_obs) - .TIE_TOL)
  list(r = r_obs, p = (cnt + 1) / (n_perm + 1), n_perm = n_perm, exact = FALSE)
}

#' Kruskal-Wallis (rank-based) multi-group association + Monte-Carlo permutation
#' @param X Numeric matrix, features x samples.
#' @param g Vector of state labels, one per sample (column of \code{X}).
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
perm_test_kruskal <- function(X, g, n_perm = 5000, seed = 42) {
  g <- factor(g); N <- length(g)
  R <- t(apply(as.matrix(X), 1, rank, ties.method = "average"))
  H_of <- function(gg) {
    M <- stats::model.matrix(~ 0 + gg)             # samples x groups
    ng <- colSums(M)
    Rbar <- sweep(R %*% M, 2, ng, "/")
    as.vector(12 / (N * (N + 1)) * ((Rbar - (N + 1) / 2)^2 %*% ng))
  }
  H <- H_of(g)
  set.seed(seed)
  cnt <- numeric(nrow(R))
  for (b in seq_len(n_perm)) cnt <- cnt + (H_of(sample(g)) >= H - .TIE_TOL)
  list(H = H, p = (cnt + 1) / (n_perm + 1), n_perm = n_perm, exact = FALSE)
}

#' Bounded adaptive FDR relaxation (PhenoRewire stats.adaptive_fdr_selection)
#' @param q Numeric vector of BH-adjusted q-values.
#' @param base_alpha Base FDR threshold.
#' @param adaptive Logical; allow bounded relaxation of the threshold.
#' @param min_selected Minimum number of features targeted by the relaxation.
#' @param alpha_ceiling Maximum FDR threshold allowed by the relaxation.
adaptive_fdr_selection <- function(q, base_alpha = 0.05, adaptive = TRUE,
                                   min_selected = 5, alpha_ceiling = 0.2) {
  fin <- is.finite(q)
  base_n <- sum(fin & q <= base_alpha)
  meta <- list(base_fdr_alpha = base_alpha, applied_fdr_alpha = base_alpha,
               adaptive_fdr_used = FALSE, selected_at_base_alpha = base_n,
               selected_feature_count = base_n)
  if (!adaptive || base_n >= min_selected) return(list(selected = fin & q <= base_alpha, meta = meta))
  cand <- sort(unique(round(c(base_alpha, alpha_ceiling, 0.01, 0.02, 0.05, 0.1, 0.2), 6)))
  cand <- cand[cand >= base_alpha & cand <= alpha_ceiling]
  counts <- vapply(cand, function(a) sum(fin & q <= a), numeric(1))
  idx <- which(counts >= min_selected)
  if (length(idx)) { a <- cand[idx[1]] } else {
    a <- cand[order(-counts, cand)][1] }
  meta$applied_fdr_alpha <- a
  meta$adaptive_fdr_used <- abs(a - base_alpha) > 1e-12
  meta$selected_feature_count <- sum(fin & q <= a)
  list(selected = fin & q <= a, meta = meta)
}

#' Feature selection for one comparison
#'
#' @param X features x samples (preprocessed)
#' @param target binary (0 = reference, 1 = case), numeric (temporal) or factor
#'   (>2 states, Kruskal-Wallis).
#' @return data.frame with r/H, p, q, selected, group means and log2FC.
#' @param mode Selection mode: "phenotype", "temporal" or "multigroup".
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param fdr_alpha Base FDR threshold for selection.
#' @param adaptive Logical; allow bounded relaxation of the threshold.
#' @param min_selected Minimum number of features targeted by the relaxation.
#' @param alpha_ceiling Maximum FDR threshold allowed by the relaxation.
#' @param seed Random seed (for reproducibility).
#' @param min_features_hard_stop Warn when fewer than this many features are selected.
select_features <- function(X, target, mode = c("phenotype", "temporal", "multigroup"),
                            n_perm = "exact", fdr_alpha = 0.05, adaptive = TRUE,
                            min_selected = 5, alpha_ceiling = 0.2, seed = 42,
                            min_features_hard_stop = 2) {
  mode <- match.arg(mode)
  res <- switch(mode,
    phenotype = perm_test_binary(X, target, n_perm = n_perm, seed = seed),
    temporal = perm_test_numeric(X, target, n_perm = ifelse(n_perm == "exact", 1000, n_perm), seed = seed),
    multigroup = perm_test_kruskal(X, target, n_perm = ifelse(n_perm == "exact", 5000, n_perm), seed = seed))
  q <- stats::p.adjust(res$p, "BH")
  sel <- adaptive_fdr_selection(q, fdr_alpha, adaptive, min_selected, alpha_ceiling)
  out <- data.frame(feature_id = rownames(X), stat = if (!is.null(res$r)) res$r else res$H,
                    p_empirical = res$p, q = q, selected = sel$selected,
                    stringsAsFactors = FALSE)
  if (mode == "phenotype") {
    out$mean_ref <- rowMeans(X[, target == 0, drop = FALSE])
    out$mean_case <- rowMeans(X[, target == 1, drop = FALSE])
    out$log2FC_case_vs_ref <- out$mean_case - out$mean_ref
  }
  meta <- c(sel$meta, list(n_features_total = nrow(X), n_perm = res$n_perm,
                           exact_permutation = res$exact, mode = mode))
  if (sum(out$selected) < min_features_hard_stop)
    warning(sprintf("Only %d features selected (hard stop %d).", sum(out$selected), min_features_hard_stop))
  attr(out, "meta") <- meta
  out
}
