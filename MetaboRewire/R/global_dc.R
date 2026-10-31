# ---- MetaboRewire - Covariate adjustment, limma selection, whole-metabolome DC --

#' Remove additive nuisance effects (e.g. processing block) while protecting the
#' biological states (equivalent to limma::removeBatchEffect with a state design).
#' @param X Numeric matrix, features x samples.
#' @param covariate Nuisance covariate (e.g. processing block), one value per sample.
#' @param state State label per sample (state effects are protected).
adjust_covariate <- function(X, covariate, state) {
  covariate <- factor(covariate); state <- factor(state)
  if (nlevels(covariate) < 2) return(X)
  D <- stats::model.matrix(~ state + covariate)
  B <- solve(crossprod(D), crossprod(D, t(X)))         # coefficients x features
  cv <- grep("^covariate", colnames(D))
  X - t(D[, cv, drop = FALSE] %*% B[cv, , drop = FALSE])
}

#' limma moderated-t selection for one two-state comparison (optional covariate)
#' @param X Numeric matrix, features x samples.
#' @param y Binary group indicator (0 = reference, 1 = case), one per sample.
#' @param covariate Nuisance covariate (e.g. processing block), one value per sample.
limma_select <- function(X, y, covariate = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("limma required")
  grp <- factor(ifelse(y == 1, "case", "ref"), levels = c("ref", "case"))
  D <- if (is.null(covariate) || length(unique(covariate)) < 2) stats::model.matrix(~ grp) else
    stats::model.matrix(~ grp + factor(covariate))
  fit <- limma::eBayes(limma::lmFit(X, D))
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  list(stat = tt$t, p = tt$P.Value, logFC = tt$logFC)
}

.rank_cor <- function(X) { Rn <- rank_rows(X); tcrossprod(Rn) }

.dc_stats <- function(Za, Zb, C) {
  Sa <- .rank_cor(Za); Sb <- .rank_cor(Zb)
  Ds <- Sb - Sa; Ra <- abs(Sa); Rb <- abs(Sb)
  diag(Ra) <- 0; diag(Rb) <- 0; diag(Ds) <- 0
  P <- nrow(Ra)
  node_dc <- rowSums(abs(Ds)) / (P - 1)
  node_dcoh <- (rowSums(Rb) - rowSums(Ra)) / (P - 1)
  # class-pair mean |rho| (pairs i != j)
  nk <- colSums(C)
  n_pair <- outer(nk, nk) - diag(nk, nrow = length(nk))  # ordered pair counts (no self)
  coh_a <- (t(C) %*% Ra %*% C) / pmax(n_pair, 1)
  coh_b <- (t(C) %*% Rb %*% C) / pmax(n_pair, 1)
  ut <- upper.tri(coh_a, diag = TRUE)
  list(node_dc = node_dc, node_dcoh = node_dcoh, coh_a = coh_a[ut], coh_b = coh_b[ut],
       global_dc = mean(abs(Ds[upper.tri(Ds)])),
       global_dcoh = mean(Rb[upper.tri(Rb)]) - mean(Ra[upper.tri(Ra)]))
}

#' Whole-metabolome, threshold-free differential co-variation between two states
#'
#' Requires no feature pre-selection and therefore also detects metabolites whose
#' mean level is unchanged. Node statistics: DC_i = mean_j |rho_B,ij - rho_A,ij|
#' and coherence change dcoh_i = mean_j (|rho_B,ij| - |rho_A,ij|). Class-pair
#' statistics: change in mean |rho| between all metabolite pairs of two pathway
#' classes (pathway cross-talk). All calibrated by exact label permutation of
#' within-state standardised data.
#' @param X Numeric matrix, features x samples.
#' @param g Vector of state labels, one per sample (column of \code{X}).
#' @param levels Character vector \code{c(reference, case)} giving the two states to compare.
#' @param classes Named character vector mapping feature id to pathway class.
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
global_dc_test <- function(X, g, levels, classes, n_perm = "exact", seed = 42) {
  g <- as.character(g)
  keep <- g %in% levels
  X <- X[, keep, drop = FALSE]; g <- g[keep]
  Z <- standardize_within(X, g)
  y <- g == levels[2]
  cl <- factor(classes[rownames(X)])
  C <- stats::model.matrix(~ 0 + cl); colnames(C) <- levels(cl)
  obs <- .dc_stats(Z[, !y], Z[, y], C)
  perms <- label_permutations(y, n_perm, seed = seed)
  exact <- isTRUE(attr(perms, "exact"))
  K <- ncol(C); ut <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  NDC <- matrix(0, nrow(perms), nrow(X)); NDH <- NDC
  NCP <- matrix(0, nrow(perms), nrow(ut)); NG <- matrix(0, nrow(perms), 2)
  for (i in seq_len(nrow(perms))) {
    s <- .dc_stats(Z[, !perms[i, ]], Z[, perms[i, ]], C)
    NDC[i, ] <- s$node_dc; NDH[i, ] <- s$node_dcoh
    NCP[i, ] <- s$coh_b - s$coh_a; NG[i, ] <- c(s$global_dc, s$global_dcoh)
  }
  if (!exact) { NDC <- rbind(obs$node_dc, NDC); NDH <- rbind(obs$node_dcoh, NDH)
    NCP <- rbind(obs$coh_b - obs$coh_a, NCP); NG <- rbind(c(obs$global_dc, obs$global_dcoh), NG) }
  p2 <- function(Nm, o) { mu <- colMeans(Nm)
    vapply(seq_along(o), function(j) mean(abs(Nm[, j] - mu[j]) >= abs(o[j] - mu[j]) - 1e-12), 1) }
  nodes <- data.frame(feature_id = rownames(X), DC = obs$node_dc, DC_null_mean = colMeans(NDC),
                      p_DC = colMeans(sweep(NDC, 2, obs$node_dc - 1e-12, ">=")),
                      dcoh = obs$node_dcoh, p_dcoh = p2(NDH, obs$node_dcoh), stringsAsFactors = FALSE)
  nodes$q_DC <- stats::p.adjust(nodes$p_DC, "BH"); nodes$q_dcoh <- stats::p.adjust(nodes$p_dcoh, "BH")
  dcp <- obs$coh_b - obs$coh_a
  nk <- colSums(C)
  npairs <- ifelse(ut[, 1] == ut[, 2], nk[ut[, 1]] * (nk[ut[, 1]] - 1) / 2, nk[ut[, 1]] * nk[ut[, 2]])
  pairs <- data.frame(class_a = colnames(C)[ut[, 1]], class_b = colnames(C)[ut[, 2]],
                      n_metabolite_pairs = npairs, coherence_a = obs$coh_a, coherence_b = obs$coh_b,
                      delta_coherence = dcp, null_sd = apply(NCP, 2, stats::sd),
                      p = p2(NCP, dcp), stringsAsFactors = FALSE)
  pairs$z <- (pairs$delta_coherence - colMeans(NCP)) / pmax(pairs$null_sd, 1e-12)
  pairs$q <- stats::p.adjust(pairs$p, "BH")
  glob <- data.frame(stat = c("global_DC", "global_delta_coherence"),
                     observed = c(obs$global_dc, obs$global_dcoh), null_mean = colMeans(NG),
                     null_sd = apply(NG, 2, stats::sd),
                     p = c(mean(NG[, 1] >= obs$global_dc - 1e-12),
                           p2(NG[, 2, drop = FALSE], obs$global_dcoh)))
  list(nodes = nodes, pairs = pairs, global = glob, exact = exact, n_perm = nrow(NDC), levels = levels)
}
