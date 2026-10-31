# ---- MetaboRewire - Pathway cross-talk and cross-compartment coupling --------

#' Pathway cross-talk enrichment of a set of edges
#'
#' Counts edges between every pair of pathway classes (including within-class)
#' and compares them with a null in which class labels are permuted among the
#' nodes of the feature universe (degree sequence of the observed network kept).
#' @param edges data.frame with node_a, node_b
#' @param classes named character vector: feature_id -> class (the universe)
#' @return long data.frame: class_a, class_b, observed, expected, log2_oe, p_enrich, p_deplete
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
pathway_crosstalk <- function(edges, classes, n_perm = 5000, seed = 42) {
  lv <- sort(unique(classes))
  K <- length(lv)
  cm <- function(cl) {
    a <- factor(cl[edges$node_a], levels = lv); b <- factor(cl[edges$node_b], levels = lv)
    M <- table(a, b); M <- M + t(M); diag(M) <- diag(M) / 2
    M[upper.tri(M, diag = TRUE)]
  }
  if (!nrow(edges)) return(data.frame())
  obs <- cm(classes)
  set.seed(seed)
  nm <- names(classes)
  null <- vapply(seq_len(n_perm), function(i) cm(stats::setNames(sample(classes), nm)), numeric(length(obs)))
  idx <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  ex <- rowMeans(null)
  data.frame(class_a = lv[idx[, 1]], class_b = lv[idx[, 2]], observed = as.numeric(obs),
             expected = ex, log2_oe = log2((as.numeric(obs) + 0.5) / (ex + 0.5)),
             p_enrich = (rowSums(null >= obs) + 1) / (n_perm + 1),
             p_deplete = (rowSums(null <= obs) + 1) / (n_perm + 1),
             stringsAsFactors = FALSE)
}

#' Same-metabolite coupling between two compartments measured in the same animals
#'
#' For every metabolite present in both matrices: Spearman rho between
#' compartments across animals (i) within each state and (ii) pooled after
#' removing state means (within-state coupling, n = all animals).
#' @param XA,XB features x samples matrices, columns named by animal id
#' @param state named vector animal id -> state
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
cross_compartment_coupling <- function(XA, XB, state, n_perm = 2000, seed = 42) {
  common_f <- intersect(rownames(XA), rownames(XB))
  ids <- intersect(colnames(XA), colnames(XB))
  A <- XA[common_f, ids, drop = FALSE]; B <- XB[common_f, ids, drop = FALSE]
  st <- as.character(state[ids])
  rowcor <- function(P, Q) rowSums(rank_rows(P) * rank_rows(Q))
  out <- data.frame(feature_id = common_f, stringsAsFactors = FALSE)
  for (s in unique(st)) out[[paste0("rho_", s)]] <- rowcor(A[, st == s], B[, st == s])
  # pooled within-state: residualise state means, then rank correlation
  res <- function(M) { for (s in unique(st)) M[, st == s] <- M[, st == s] - rowMeans(M[, st == s, drop = FALSE]); M }
  Ar <- res(A); Br <- res(B)
  out$rho_within_pooled <- rowcor(Ar, Br)
  # permutation p for pooled coupling: permute animals within state in one compartment
  set.seed(seed)
  cnt <- numeric(nrow(out))
  for (b in seq_len(n_perm)) {
    perm <- seq_along(ids)
    for (s in unique(st)) { w <- which(st == s); perm[w] <- w[sample.int(length(w))] }
    cnt <- cnt + (abs(rowcor(Ar, Br[, perm, drop = FALSE])) >= abs(out$rho_within_pooled) - 1e-12)
  }
  out$p_pooled <- (cnt + 1) / (n_perm + 1)
  out$q_pooled <- stats::p.adjust(out$p_pooled, "BH")
  out
}

# Bipartite (cross-compartment) edges with BH over all pairs
.bipartite_keys <- function(A, B, th) {
  n <- ncol(A)
  R <- suppressWarnings(stats::cor(t(A), t(B), method = "spearman"))
  R[!is.finite(R)] <- 0
  tt <- R * sqrt((n - 2) / pmax(1 - R^2, 0)); P <- 2 * stats::pt(-abs(tt), n - 2)
  P[abs(R) >= 1 - 1e-12] <- 0
  q <- matrix(stats::p.adjust(P, "BH"), nrow(P))
  w <- which(q <= th$fdr & abs(R) >= th$abs_r)
  list(key = w, sign = sign(R[w]), r = R[w], q = q[w])
}

#' Rewiring of the bipartite inter-compartment network between two states,
#' with a paired label-permutation null (animals keep both compartments).
#' @param A Numeric matrix, features x animals, compartment A.
#' @param B Numeric matrix, features x animals, compartment B (same animals and order as \code{A}).
#' @param g Vector of state labels, one per sample (column of \code{X}).
#' @param th List with \code{abs_r} and \code{fdr} thresholds for bipartite edges.
#' @param n_perm Number of permutations, or "exact" to enumerate all labelings when feasible.
#' @param seed Random seed (for reproducibility).
#' @param min_r Minimum absolute correlation in both states for a sign switch.
bipartite_rewiring <- function(A, B, g, th = list(abs_r = 0.8, fdr = 0.2), n_perm = "exact",
                               seed = 42, min_r = 0.5) {
  g <- as.character(g); lv <- unique(g); y <- g == lv[2]
  ZA <- standardize_within(A, g); ZB <- standardize_within(B, g)
  stat <- function(mask) {
    ka <- .bipartite_keys(ZA[, !mask], ZB[, !mask], th)
    kb <- .bipartite_keys(ZA[, mask], ZB[, mask], th)
    sh <- intersect(ka$key, kb$key)
    c(n_a = length(ka$key), n_b = length(kb$key), n_shared = length(sh),
      n_specific = length(setdiff(union(ka$key, kb$key), sh)))
  }
  obs <- stat(y)
  perms <- label_permutations(y, n_perm, seed = seed)
  null <- t(vapply(seq_len(nrow(perms)), function(i) stat(perms[i, ]), numeric(4)))
  if (!isTRUE(attr(perms, "exact"))) null <- rbind(obs, null)
  ka <- .bipartite_keys(A[, !y], B[, !y], th); kb <- .bipartite_keys(A[, y], B[, y], th)
  edge_df <- function(k, lab) {
    if (!length(k$key)) return(NULL)
    ij <- arrayInd(k$key, c(nrow(A), nrow(B)))
    data.frame(feature_A = rownames(A)[ij[, 1]], feature_B = rownames(B)[ij[, 2]],
               r = k$r, q = k$q, state = lab, stringsAsFactors = FALSE)
  }
  list(observed = obs,
       global = data.frame(stat = names(obs), observed = obs, null_mean = colMeans(null),
                           p_upper = colMeans(sweep(null, 2, obs, ">=")), stringsAsFactors = FALSE),
       edges = rbind(edge_df(ka, lv[1]), edge_df(kb, lv[2])), labels = lv)
}
