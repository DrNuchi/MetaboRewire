# ---- MetaboRewire - Three-state trajectory rewiring (extension) --------------
# Design: reference (R, e.g. healthy) -> perturbation (D, e.g. EAM) ->
# intervention (T, e.g. EAM + metformin). Each edge is assigned a presence
# pattern over (R, D, T) and a biological class:
#   conserved                 1 1 1
#   disease_gained_rescued    0 1 0   gained in disease, lost again under treatment
#   disease_gained_persistent 0 1 1
#   disease_lost_restored     1 0 1   lost in disease, restored under treatment
#   disease_lost_persistent   1 0 0
#   treatment_emergent        0 0 1
#   treatment_disrupted       1 1 0
# Edges present in >1 state with opposite signs are flagged sign_switch.
# Network rescue index (NRI) = (gained_rescued + lost_restored) /
#                              (all disease-perturbed edges).

.joint_thresholds <- function(Ulist, p) {
  G <- threshold_grid(p$abs_r, p$fdr, p$abs_r_floor, p$fdr_ceiling)
  if (!p$adaptive) return(list(abs_r = p$abs_r, fdr = p$fdr))
  sc <- vapply(seq_len(nrow(G)), function(i) {
    sum(vapply(Ulist, function(U) {
      if (is.null(U)) return(0)
      fin <- is.finite(U$r) & is.finite(U$q)
      ec <- sum(fin & U$q <= G$fdr[i] & abs(U$r) >= G$abs_r[i])
      .threshold_score(ec, U$n_nodes * (U$n_nodes - 1) / 2, G$abs_r[i], G$fdr[i], p)
    }, numeric(1)))
  }, numeric(1))
  list(abs_r = G$abs_r[which.min(sc)], fdr = G$fdr[which.min(sc)])
}

#' Classify edges by presence pattern across three state networks
#' @param keysets List of three edge sets (elements \code{key}, \code{r}), one per state.
#' @param states Character vector of state labels in trajectory order (reference, perturbation, intervention).
#' @param min_r Minimum absolute correlation in both states for a sign switch.
classify_three_state_edges <- function(keysets, states = c("R", "D", "T"), min_r = 0.5) {
  all_k <- sort(unique(unlist(lapply(keysets, `[[`, "key"))))
  nk <- length(all_k)
  pres <- matrix(unlist(lapply(keysets, function(k) all_k %in% k$key)), nk, length(keysets))
  rmat <- matrix(unlist(lapply(keysets, function(k) as.numeric(k$r[match(all_k, k$key)]))), nk, length(keysets))
  pat <- apply(pres * 1L, 1, paste, collapse = "")
  cls <- c("111" = "conserved", "010" = "disease_gained_rescued", "011" = "disease_gained_persistent",
           "101" = "disease_lost_restored", "100" = "disease_lost_persistent",
           "001" = "treatment_emergent", "110" = "treatment_disrupted")[pat]
  sgn <- apply(rmat, 1, function(r) { r <- r[!is.na(r) & abs(r) >= min_r]
    length(r) > 1 && length(unique(sign(r))) > 1 })
  cls <- unname(cls); cls[sgn] <- "sign_switch"
  data.frame(key = all_k, pattern = pat, class = cls, sign_switch = sgn,
             r_R = rmat[, 1], r_D = rmat[, 2], r_T = rmat[, 3], stringsAsFactors = FALSE)
}

.class_counts <- function(ed) {
  lv <- c("conserved", "disease_gained_rescued", "disease_gained_persistent",
          "disease_lost_restored", "disease_lost_persistent", "treatment_emergent",
          "treatment_disrupted", "sign_switch")
  cc <- table(factor(ed$class, levels = lv))
  pert <- sum(cc[c("disease_gained_rescued", "disease_gained_persistent",
                   "disease_lost_restored", "disease_lost_persistent")])
  c(cc,
    disease_perturbed = pert,
    NRI = if (pert > 0) unname((cc["disease_gained_rescued"] + cc["disease_lost_restored"]) / pert) else NA,
    treatment_specific = unname(cc["treatment_emergent"] + cc["treatment_disrupted"]))
}

#' Three-state rewiring with jointly chosen thresholds
#'
#' @param X features x samples
#' @param g state labels; `states` gives (reference, perturbation, intervention)
#' @param states Character vector of state labels in trajectory order (reference, perturbation, intervention).
#' @param p List of network parameters from \code{network_params()}.
#' @param feature_names Optional feature names (currently unused).
three_state_rewiring <- function(X, g, states, p = network_params(), feature_names = NULL) {
  g <- as.character(g)
  Ul <- lapply(states, function(s) .upper_edges(X[, g == s, drop = FALSE]))
  th <- .joint_thresholds(Ul, p)
  ks <- lapply(states, function(s) .net_keys(X[, g == s, drop = FALSE], p, fixed = th))
  ed <- classify_three_state_edges(ks, states, p$sign_switch_min_r)
  P <- nrow(X)
  ed$node_a <- rownames(X)[(ed$key - 1L) %/% P + 1L]
  ed$node_b <- rownames(X)[(ed$key - 1L) %% P + 1L]
  counts <- .class_counts(ed)
  # node-level rescue profile
  long <- rbind(data.frame(f = ed$node_a, c = ed$class), data.frame(f = ed$node_b, c = ed$class))
  tab <- as.data.frame.matrix(table(long$f, factor(long$c, levels = names(counts)[1:8])))
  tab$feature_id <- rownames(tab)
  tab$disease_perturbed <- rowSums(tab[, c("disease_gained_rescued", "disease_gained_persistent",
                                           "disease_lost_restored", "disease_lost_persistent")])
  tab$node_rescue_fraction <- ifelse(tab$disease_perturbed > 0,
    (tab$disease_gained_rescued + tab$disease_lost_restored) / tab$disease_perturbed, NA)
  list(edges = ed, counts = counts, thresholds = th, nodes = tab, states = states,
       n_edges_state = sapply(ks, function(k) length(k$key)))
}

#' Null distribution of the network rescue index
#'
#' H0: the intervention does not change the dependence structure relative to the
#' perturbation state, i.e. D and T samples are exchangeable after within-state
#' standardisation. The reference network is fixed; D/T labels are permuted
#' (exactly when feasible) and the classification recomputed with the observed
#' thresholds.
#' @param X Numeric matrix, features x samples.
#' @param g Vector of state labels, one per sample (column of \code{X}).
#' @param states Character vector of state labels in trajectory order (reference, perturbation, intervention).
#' @param obs Result of \code{three_state_rewiring()} on the observed data.
#' @param p List of network parameters from \code{network_params()}.
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
rescue_permutation_test <- function(X, g, states, obs, p = network_params(),
                                    n_perm = "exact", seed = 42) {
  g <- as.character(g)
  sub <- g %in% states[2:3]
  Z <- X
  Z[, sub] <- standardize_within(X[, sub, drop = FALSE], g[sub])
  kR <- .net_keys(X[, g == states[1], drop = FALSE], p, fixed = obs$thresholds)
  Zs <- Z[, sub, drop = FALSE]; ys <- g[sub] == states[3]
  perms <- label_permutations(ys, n_perm, seed = seed)
  null <- t(vapply(seq_len(nrow(perms)), function(i) {
    m <- perms[i, ]
    kD <- .net_keys(Zs[, !m, drop = FALSE], p, fixed = obs$thresholds)
    kT <- .net_keys(Zs[, m, drop = FALSE], p, fixed = obs$thresholds)
    .class_counts(classify_three_state_edges(list(kR, kD, kT), states, p$sign_switch_min_r))
  }, numeric(length(obs$counts))))
  if (!isTRUE(attr(perms, "exact"))) null <- rbind(obs$counts, null)
  stats <- c("disease_gained_rescued", "disease_lost_restored", "treatment_emergent",
             "treatment_disrupted", "treatment_specific", "NRI", "conserved")
  out <- data.frame(stat = stats, observed = as.numeric(obs$counts[stats]),
                    null_mean = colMeans(null[, stats], na.rm = TRUE),
                    null_sd = apply(null[, stats], 2, stats::sd, na.rm = TRUE),
                    p_upper = vapply(stats, function(s) mean(null[, s] >= obs$counts[s], na.rm = TRUE), 1),
                    p_lower = vapply(stats, function(s) mean(null[, s] <= obs$counts[s], na.rm = TRUE), 1),
                    stringsAsFactors = FALSE)
  out$z <- (out$observed - out$null_mean) / pmax(out$null_sd, 1e-12)
  list(table = out, null = null, exact = isTRUE(attr(perms, "exact")))
}

#' Abundance-level rescue for features changed between reference and perturbation
#'
#' rescue = 1 - (mean_T - mean_R) / (mean_D - mean_R)  (log2 scale);
#' 1 = full normalisation, 0 = no effect, < 0 = exacerbation, > 1 = overshoot.
#' @param X Numeric matrix, features x samples.
#' @param g Vector of state labels, one per sample (column of \code{X}).
#' @param states Character vector of state labels in trajectory order (reference, perturbation, intervention).
abundance_rescue <- function(X, g, states) {
  m <- sapply(states, function(s) rowMeans(X[, g == s, drop = FALSE]))
  d <- m[, 2] - m[, 1]; t <- m[, 3] - m[, 1]
  data.frame(feature_id = rownames(X), log2FC_D_vs_R = d, log2FC_T_vs_R = t,
             log2FC_T_vs_D = m[, 3] - m[, 2], abundance_rescue = 1 - t / d,
             stringsAsFactors = FALSE)
}

#' Threshold-free three-state rescue analysis (network and abundance level)
#'
#' Node coherence c_s(i) = mean_j |rho_s(i, j)| is computed in every state.
#' Disease effect d1 = c_D - c_R, treatment effect d2 = c_T - c_D. A treatment
#' that reverses the disease-associated change produces a negative association
#' between d1 and d2 across metabolites. Because c_D enters both differences,
#' sampling noise alone induces a negative association (regression to the
#' mean); the statistic is therefore calibrated against the null in which D and
#' T samples are exchangeable (labels permuted after within-state
#' standardisation for the network level, on the raw data for the abundance
#' level). Rescue indices are reported as the observed association minus its
#' null expectation.
#' @param classes named vector feature -> pathway class (for class-pair rescue)
#' @param X Numeric matrix, features x samples.
#' @param g Vector of state labels, one per sample (column of \code{X}).
#' @param states Character vector of state labels in trajectory order (reference, perturbation, intervention).
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
three_state_dc <- function(X, g, states, classes, n_perm = "exact", seed = 42) {
  g <- as.character(g)
  keep <- g %in% states; X <- X[, keep, drop = FALSE]; g <- g[keep]
  cl <- factor(classes[rownames(X)]); C <- stats::model.matrix(~ 0 + cl); colnames(C) <- levels(cl)
  nk <- colSums(C); n_pair <- outer(nk, nk) - diag(nk, nrow = length(nk))
  K <- ncol(C); ut <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  coh <- function(Z) { S <- abs(tcrossprod(rank_rows(Z))); diag(S) <- 0
    list(node = rowSums(S) / (nrow(S) - 1), pair = ((t(C) %*% S %*% C) / pmax(n_pair, 1))[ut]) }
  cR <- coh(X[, g == states[1]])
  sub <- g %in% states[2:3]
  Zs <- standardize_within(X[, sub], g[sub]); Xs <- X[, sub]; ys <- g[sub] == states[3]
  mR <- rowMeans(X[, g == states[1], drop = FALSE])
  stat <- function(mask) {
    cD <- coh(Zs[, !mask]); cT <- coh(Zs[, mask])
    d1 <- cD$node - cR$node; d2 <- cT$node - cD$node
    p1 <- cD$pair - cR$pair; p2 <- cT$pair - cD$pair
    a1 <- rowMeans(Xs[, !mask, drop = FALSE]) - mR; a2 <- rowMeans(Xs[, mask, drop = FALSE]) - rowMeans(Xs[, !mask, drop = FALSE])
    list(net_cor = stats::cor(d1, d2), abund_cor = stats::cor(a1, a2),
         pair_d1 = p1, pair_d2 = p2, node_d1 = d1, node_d2 = d2, cD = cD, cT = cT)
  }
  obs <- stat(ys)
  perms <- label_permutations(ys, n_perm, seed = seed)
  nl <- t(vapply(seq_len(nrow(perms)), function(i) {
    s <- stat(perms[i, ]); c(s$net_cor, s$abund_cor, s$pair_d1 * s$pair_d2) }, numeric(2 + nrow(ut))))
  if (!isTRUE(attr(perms, "exact"))) nl <- rbind(c(obs$net_cor, obs$abund_cor, obs$pair_d1 * obs$pair_d2), nl)
  glob <- data.frame(level = c("network_coherence", "abundance"),
                     observed_cor = c(obs$net_cor, obs$abund_cor),
                     null_mean_cor = colMeans(nl[, 1:2]),
                     rescue_index = -(c(obs$net_cor, obs$abund_cor) - colMeans(nl[, 1:2])),
                     p_rescue = c(mean(nl[, 1] <= obs$net_cor + 1e-12), mean(nl[, 2] <= obs$abund_cor + 1e-12)))
  prod_obs <- obs$pair_d1 * obs$pair_d2
  pairs <- data.frame(class_a = colnames(C)[ut[, 1]], class_b = colnames(C)[ut[, 2]],
                      coherence_R = cR$pair, coherence_D = obs$cD$pair, coherence_T = obs$cT$pair,
                      delta_disease = obs$pair_d1, delta_treatment = obs$pair_d2,
                      reversal_product = prod_obs,
                      p_reversal = vapply(seq_along(prod_obs), function(j) mean(nl[, 2 + j] <= prod_obs[j] + 1e-12), 1),
                      stringsAsFactors = FALSE)
  pairs$pattern <- ifelse(pairs$delta_disease < 0 & pairs$delta_treatment > 0, "decoupled_then_recoupled",
                   ifelse(pairs$delta_disease > 0 & pairs$delta_treatment < 0, "hypercoupled_then_normalised",
                   ifelse(pairs$delta_disease < 0, "decoupled_persistent_or_worse", "hypercoupled_persistent_or_worse")))
  pairs$q_reversal <- stats::p.adjust(pairs$p_reversal, "BH")
  nodes <- data.frame(feature_id = rownames(X), coherence_R = cR$node, coherence_D = obs$cD$node,
                      coherence_T = obs$cT$node, delta_disease = obs$node_d1,
                      delta_treatment = obs$node_d2, stringsAsFactors = FALSE)
  list(global = glob, pairs = pairs, nodes = nodes, null = nl[, 1:2], exact = isTRUE(attr(perms, "exact")))
}
