# ---- MetaboRewire - State-specific co-variation networks ---------------------
# Port of PhenoRewire networks.py: Spearman correlation network per state with
# BH-adjusted edge p-values, |rho| threshold and bounded adaptive thresholding
# (target: >= min_edges and <= max_density), plus igraph topology metrics.

#' Spearman correlation matrix with t-approximation p-values (as scipy.stats.spearmanr)
#' @param X Numeric matrix, features x samples.
spearman_cor_p <- function(X) {
  n <- ncol(X)
  R <- suppressWarnings(stats::cor(t(X), method = "spearman"))
  tt <- R * sqrt((n - 2) / pmax(1 - R^2, 0))
  P <- 2 * stats::pt(-abs(tt), df = n - 2)
  P[abs(R) >= 1 - 1e-12] <- 0
  list(R = R, P = P, n = n)
}

#' Adaptive threshold grid (PhenoRewire networks._threshold_grid)
#' @param base_abs_r Base absolute correlation threshold.
#' @param base_fdr Base FDR threshold.
#' @param floor_abs_r Lowest absolute correlation threshold in the grid.
#' @param ceil_fdr Highest FDR threshold in the grid.
threshold_grid <- function(base_abs_r, base_fdr, floor_abs_r, ceil_fdr) {
  a <- sort(unique(round(c(min(0.95, base_abs_r + 0.10), min(0.95, base_abs_r + 0.05),
                           base_abs_r, max(floor_abs_r, base_abs_r - 0.05),
                           max(floor_abs_r, base_abs_r - 0.10), floor_abs_r), 4)),
            decreasing = TRUE)
  f <- sort(unique(round(c(min(base_fdr, 0.01), base_fdr,
                           min(ceil_fdr, max(base_fdr, 0.10)), ceil_fdr), 4)))
  expand <- expand.grid(fdr = f, abs_r = a)[, c("abs_r", "fdr")]
  expand
}

.threshold_score <- function(edge_count, total_possible, abs_thr, fdr_thr, p) {
  tmin <- min(p$min_edges, total_possible)
  tmax <- max(tmin, floor(p$max_density * total_possible))
  too_few <- pmax(0, tmin - edge_count); too_many <- pmax(0, edge_count - tmax)
  dev <- abs(abs_thr - p$abs_r) / 0.05 + abs(fdr_thr - p$fdr) / max(p$fdr, 1e-6)
  too_few * 1000 + too_many * 10 + dev + edge_count / total_possible
}

#' Choose (|rho|, FDR) thresholds (PhenoRewire networks._select_adaptive_thresholds)
#' @param r Numeric vector of edge correlations.
#' @param q Numeric vector of BH-adjusted q-values.
#' @param n_nodes Number of network nodes.
#' @param p List of network parameters from \code{network_params()}.
select_adaptive_thresholds <- function(r, q, n_nodes, p) {
  total <- n_nodes * (n_nodes - 1) / 2
  fin <- is.finite(r) & is.finite(q)
  if (!p$adaptive || total == 0 || sum(fin) == 0)
    return(list(abs_r = p$abs_r, fdr = p$fdr, adaptive_used = FALSE))
  G <- threshold_grid(p$abs_r, p$fdr, p$abs_r_floor, p$fdr_ceiling)
  best <- NULL; best_score <- Inf
  for (i in seq_len(nrow(G))) {
    ec <- sum(fin & q <= G$fdr[i] & abs(r) >= G$abs_r[i])
    s <- .threshold_score(ec, total, G$abs_r[i], G$fdr[i], p)
    if (s < best_score) { best_score <- s; best <- i }
  }
  list(abs_r = G$abs_r[best], fdr = G$fdr[best],
       adaptive_used = abs(G$abs_r[best] - p$abs_r) > 1e-12 || abs(G$fdr[best] - p$fdr) > 1e-12)
}

#' Default network parameters (PhenoRewire phenotype defaults used in the paper)
#' @param abs_r Absolute Spearman correlation threshold.
#' @param fdr Edge-level BH FDR threshold.
#' @param adaptive Logical; allow bounded relaxation of the threshold.
#' @param min_edges Target minimum number of edges for adaptive thresholds.
#' @param max_density Target maximum network density for adaptive thresholds.
#' @param abs_r_floor Lowest absolute correlation threshold allowed by adaptive selection.
#' @param fdr_ceiling Highest FDR threshold allowed by adaptive selection.
#' @param sign_switch_min_r Minimum absolute correlation in both networks for an edge to count as a sign switch.
network_params <- function(abs_r = 0.70, fdr = 0.05, adaptive = TRUE, min_edges = 3,
                           max_density = 0.25, abs_r_floor = 0.50, fdr_ceiling = 0.20,
                           sign_switch_min_r = 0.5) {
  list(abs_r = abs_r, fdr = fdr, adaptive = adaptive, min_edges = min_edges,
       max_density = max_density, abs_r_floor = abs_r_floor, fdr_ceiling = fdr_ceiling,
       sign_switch_min_r = sign_switch_min_r)
}

.upper_edges <- function(X) {
  sd_ok <- apply(X, 1, stats::sd) > 0
  X2 <- X[sd_ok, , drop = FALSE]
  if (nrow(X2) < 2 || ncol(X2) < 3) return(NULL)
  cp <- spearman_cor_p(X2)
  iu <- which(upper.tri(cp$R), arr.ind = TRUE)
  r <- cp$R[iu]; pv <- cp$P[iu]
  fin <- is.finite(r) & is.finite(pv)
  q <- rep(NA_real_, length(pv)); q[fin] <- stats::p.adjust(pv[fin], "BH")
  list(feats = rownames(X2), i = iu[, 1], j = iu[, 2], r = r, q = q, n_nodes = nrow(X2))
}

#' Build one state-specific correlation network
#'
#' @param X features x samples (only samples of one state)
#' @param p list from network_params()
#' @param fixed optional list(abs_r, fdr) to bypass adaptive selection
#' @return list(edges = data.frame(node_a, node_b, r, q, sign), thresholds, n_nodes)
correlation_network <- function(X, p = network_params(), fixed = NULL) {
  empty <- data.frame(node_a = character(), node_b = character(), r = numeric(),
                      q = numeric(), sign = character(), stringsAsFactors = FALSE)
  U <- .upper_edges(X)
  if (is.null(U)) return(list(edges = empty, thresholds = list(abs_r = p$abs_r, fdr = p$fdr,
                                                               adaptive_used = FALSE), n_nodes = 0))
  th <- if (is.null(fixed)) select_adaptive_thresholds(U$r, U$q, U$n_nodes, p) else
    c(fixed, list(adaptive_used = NA))
  keep <- is.finite(U$r) & is.finite(U$q) & U$q <= th$fdr & abs(U$r) >= th$abs_r
  e <- data.frame(node_a = U$feats[U$i[keep]], node_b = U$feats[U$j[keep]],
                  r = U$r[keep], q = U$q[keep], stringsAsFactors = FALSE)
  e$sign <- ifelse(e$r >= 0, "pos", "neg")
  th$n_edges <- nrow(e)
  th$density <- if (U$n_nodes > 1) nrow(e) / (U$n_nodes * (U$n_nodes - 1) / 2) else 0
  list(edges = if (nrow(e)) e else empty, thresholds = th, n_nodes = U$n_nodes)
}

#' Topological metrics for a state network (degree, centralities, Louvain, Q)
#' @param nodes Character vector of all node (feature) ids.
#' @param edges Edge data.frame (node_a, node_b, r) of the observed network.
#' @param seed Random seed (for reproducibility).
network_metrics <- function(nodes, edges, seed = 42) {
  g <- igraph::graph_from_data_frame(
    if (nrow(edges)) data.frame(from = edges$node_a, to = edges$node_b, r = edges$r) else
      data.frame(from = character(), to = character(), r = numeric()),
    directed = FALSE, vertices = data.frame(name = unique(nodes)))
  w <- abs(igraph::E(g)$r)
  n <- igraph::vcount(g)
  deg <- igraph::degree(g)
  out <- data.frame(feature_id = igraph::V(g)$name, degree = as.numeric(deg),
                    degree_centrality = if (n > 1) deg / (n - 1) else 0,
                    weighted_degree_abs = if (length(w)) igraph::strength(g, weights = w) else 0,
                    stringsAsFactors = FALSE)
  if (igraph::ecount(g) > 0) {
    out$betweenness_centrality <- igraph::betweenness(g, weights = 1 / pmax(w, 1e-12), normalized = TRUE)
    ev <- rep(0, n); comp <- igraph::components(g)
    for (cid in unique(comp$membership)) {
      idx <- which(comp$membership == cid)
      if (length(idx) < 2) next
      sg <- igraph::induced_subgraph(g, idx)
      v <- igraph::eigen_centrality(sg, weights = abs(igraph::E(sg)$r))$vector
      ev[idx] <- v / sqrt(sum(v^2))      # unit L2 norm per component, as networkx
    }
    out$eigenvector_centrality <- ev
    set.seed(seed)
    cl <- igraph::cluster_louvain(g, weights = w)
    memb <- igraph::membership(cl)
    Q <- igraph::modularity(g, memb, weights = w)
  } else {
    out$betweenness_centrality <- 0; out$eigenvector_centrality <- 0
    memb <- seq_len(n); Q <- 0
  }
  out$community_id <- as.integer(memb)
  out$community_size <- as.integer(table(memb)[as.character(memb)])
  comp <- igraph::components(g)
  summary <- data.frame(n_nodes = n, n_edges = igraph::ecount(g),
                        density = if (n > 1) igraph::edge_density(g) else 0,
                        n_components = comp$no, largest_component_size = max(comp$csize, 0),
                        positive_edges = sum(edges$sign == "pos"), negative_edges = sum(edges$sign == "neg"),
                        mean_abs_r = if (nrow(edges)) mean(abs(edges$r)) else 0,
                        n_communities = length(unique(memb[deg > 0])),
                        modularity = Q, low_modularity_flag = Q < 0.3)
  list(nodes = out, summary = summary, graph = g)
}
