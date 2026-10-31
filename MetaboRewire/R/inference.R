# ---- MetaboRewire - Permutation inference for rewiring (extension) -----------
# PhenoRewire reports rewiring counts and scores without a null model. With
# n = 4-6 replicates per state, state-specific Spearman networks are noisy and
# two independently estimated networks differ substantially even when the true
# correlation structure is identical. MetaboRewire therefore calibrates every
# rewiring statistic against a label-permutation null:
#   1. each feature is standardised within its true state (removes location and
#      scale differences, i.e. the differential-abundance signal used for
#      feature selection, leaving only the dependence structure);
#   2. state labels are permuted (all choose(n, n_a) assignments when feasible);
#   3. the *entire* network procedure, including adaptive thresholding, is
#      re-run on each pseudo-partition;
#   4. empirical p-values are computed for global and node-level statistics.

#' Standardise each feature within each level of `g` (mean 0, sd 1)
#' @param X Numeric matrix, features x samples.
#' @param g Vector of state labels, one per sample (column of \code{X}).
standardize_within <- function(X, g) {
  Z <- X
  for (lv in unique(g)) {
    idx <- which(g == lv)
    m <- rowMeans(X[, idx, drop = FALSE])
    s <- apply(X[, idx, drop = FALSE], 1, stats::sd); s[!is.finite(s) | s == 0] <- 1
    Z[, idx] <- (X[, idx, drop = FALSE] - m) / s
  }
  Z
}

# Edge keys (integer, over the global feature index) of one network
.net_keys <- function(X, p, fixed = NULL) {
  U <- .upper_edges(X)
  if (is.null(U)) return(list(key = integer(), sign = integer(), r = numeric(), th = NULL))
  th <- if (is.null(fixed)) select_adaptive_thresholds(U$r, U$q, U$n_nodes, p) else fixed
  keep <- is.finite(U$r) & is.finite(U$q) & U$q <= th$fdr & abs(U$r) >= th$abs_r
  gi <- match(U$feats, rownames(X))
  a <- gi[U$i[keep]]; b <- gi[U$j[keep]]
  P <- nrow(X)
  list(key = (pmin(a, b) - 1L) * P + pmax(a, b), sign = sign(U$r[keep]), r = U$r[keep], th = th)
}

# Global and node statistics from two key sets
.rewire_stats <- function(A, B, P, min_r) {
  sh <- intersect(A$key, B$key)
  ao <- setdiff(A$key, B$key); bo <- setdiff(B$key, A$key)
  ia <- match(sh, A$key); ib <- match(sh, B$key)
  ss <- sh[A$sign[ia] != B$sign[ib] & abs(A$r[ia]) >= min_r & abs(B$r[ib]) >= min_r]
  ends <- function(k) c((k - 1L) %/% P + 1L, (k - 1L) %% P + 1L)
  node <- tabulate(ends(c(ao, bo, ss)), nbins = P)
  c(n_a = length(A$key), n_b = length(B$key), n_shared = length(sh),
    n_specific = length(ao) + length(bo), n_sign_switch = length(ss),
    rewiring_proportion = (length(ao) + length(bo)) / max(1, length(union(A$key, B$key))),
    jaccard_distance = 1 - length(sh) / max(1, length(union(A$key, B$key))),
    node)
}

# Threshold-free differential co-variation (MetaboRewire): uses the full
# Spearman matrices rather than thresholded edges.
.soft_stats <- function(Xa, Xb) {
  Ra <- suppressWarnings(stats::cor(t(Xa), method = "spearman")); Ra[!is.finite(Ra)] <- 0
  Rb <- suppressWarnings(stats::cor(t(Xb), method = "spearman")); Rb[!is.finite(Rb)] <- 0
  ut <- upper.tri(Ra)
  D <- abs(Ra - Rb); diag(D) <- 0
  c(DC_global = mean(D[ut]), coherence_a = mean(abs(Ra[ut])), coherence_b = mean(abs(Rb[ut])),
    delta_coherence = mean(abs(Rb[ut])) - mean(abs(Ra[ut])),
    rowSums(D) / (nrow(D) - 1))
}

#' Label-permutation test for two-state rewiring
#'
#' Hard (PhenoRewire edge classification) and soft (threshold-free differential
#' connectivity, DC) statistics are computed on the observed partition and on
#' every permuted partition of within-state standardised data.
#'
#' @param X features x samples (selected features; samples of both states)
#' @param g two-level vector of state labels
#' @param p network_params()
#' @param n_perm "exact" or integer
#' @param levels c(reference, case); default order of appearance
#' @return list(global, nodes, null, exact, n_perm, labels)
#' @param seed Random seed (for reproducibility).
#' @param max_exact Maximum number of labelings for which exact enumeration is used.
#' @param verbose Logical; print progress messages.
rewiring_permutation_test <- function(X, g, p = network_params(), n_perm = "exact",
                                      seed = 42, max_exact = 20000, verbose = FALSE,
                                      levels = NULL) {
  g <- as.character(g); lv <- if (is.null(levels)) unique(g) else levels
  stopifnot(length(lv) == 2)
  Z <- standardize_within(X, g)
  y <- g == lv[2]
  P <- nrow(X)
  hard_names <- c("n_a", "n_b", "n_shared", "n_specific", "n_sign_switch",
                  "rewiring_proportion", "jaccard_distance")
  soft_names <- c("DC_global", "coherence_a", "coherence_b", "delta_coherence")
  run <- function(mask) {
    Xa <- Z[, !mask, drop = FALSE]; Xb <- Z[, mask, drop = FALSE]
    h <- .rewire_stats(.net_keys(Xa, p), .net_keys(Xb, p), P, p$sign_switch_min_r)
    s <- .soft_stats(Xa, Xb)
    c(h[hard_names], s[soft_names], h[-(seq_along(hard_names))], s[-(seq_along(soft_names))])
  }
  obs <- run(y)
  perms <- label_permutations(y, n_perm, max_exact, seed)
  exact <- isTRUE(attr(perms, "exact"))
  null <- t(vapply(seq_len(nrow(perms)), function(i) {
    if (verbose && i %% 100 == 0) message("perm ", i, "/", nrow(perms))
    run(perms[i, ])
  }, numeric(length(obs))))
  if (!exact) null <- rbind(obs, null)
  gn <- c(hard_names, soft_names)
  G <- null[, gn, drop = FALSE]
  d_obs <- obs["n_b"] - obs["n_a"]; d_null <- G[, "n_b"] - G[, "n_a"]
  two_sided <- c("delta_coherence")
  pu <- colMeans(sweep(G, 2, obs[gn], ">=") )
  pt <- vapply(gn, function(s) mean(abs(G[, s] - mean(G[, s])) >= abs(obs[s] - mean(G[, s])) - 1e-12), 1)
  glob <- data.frame(
    stat = c(gn, "delta_edges_b_minus_a"),
    observed = c(obs[gn], d_obs),
    null_mean = c(colMeans(G), mean(d_null)),
    null_sd = c(apply(G, 2, stats::sd), stats::sd(d_null)),
    p_value = c(ifelse(gn %in% two_sided, pt, pu), mean(abs(d_null) >= abs(d_obs))),
    stringsAsFactors = FALSE)
  glob$z <- (glob$observed - glob$null_mean) / pmax(glob$null_sd, 1e-12)
  glob$test <- ifelse(glob$stat %in% c(two_sided, "delta_edges_b_minus_a"), "two-sided", "upper-tail")
  k <- length(gn)
  Nh <- null[, (k + 1):(k + P), drop = FALSE]; Ns <- null[, (k + P + 1):(k + 2 * P), drop = FALSE]
  oh <- obs[(k + 1):(k + P)]; os <- obs[(k + P + 1):(k + 2 * P)]
  ph <- colMeans(sweep(Nh, 2, oh, ">=")); ps <- colMeans(sweep(Ns, 2, os - 1e-12, ">="))
  nodes <- data.frame(feature_id = rownames(X), rewiring_score = as.numeric(oh),
                      null_mean = colMeans(Nh), p_perm = ph, q_perm = stats::p.adjust(ph, "BH"),
                      DC = as.numeric(os), DC_null_mean = colMeans(Ns), p_DC = ps,
                      q_DC = stats::p.adjust(ps, "BH"), stringsAsFactors = FALSE)
  nodes <- nodes[order(nodes$p_DC, -nodes$DC), ]
  list(global = glob, nodes = nodes, null = G, exact = exact, n_perm = nrow(null), labels = lv)
}

#' Leave-one-out support of each observed edge
#'
#' Support = fraction of leave-one-sample-out replicates in which the edge keeps
#' its sign and |rho| >= the applied threshold.
#' @param X Numeric matrix, features x samples.
#' @param edges Edge data.frame (node_a, node_b, r) of the observed network.
#' @param abs_r Absolute Spearman correlation threshold.
jackknife_edge_support <- function(X, edges, abs_r) {
  if (!nrow(edges)) return(numeric())
  n <- ncol(X)
  hits <- numeric(nrow(edges))
  ia <- match(edges$node_a, rownames(X)); ib <- match(edges$node_b, rownames(X))
  for (k in seq_len(n)) {
    R <- suppressWarnings(stats::cor(t(X[, -k, drop = FALSE]), method = "spearman"))
    rr <- R[cbind(ia, ib)]
    hits <- hits + (is.finite(rr) & abs(rr) >= abs_r & sign(rr) == sign(edges$r))
  }
  hits / n
}
