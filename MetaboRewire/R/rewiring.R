# ---- MetaboRewire - Edge-level rewiring and node triage ----------------------
# Port of PhenoRewire compute_rewiring / summarize_rewiring / triage.

.edge_key <- function(a, b) ifelse(a <= b, paste(a, b, sep = "__"), paste(b, a, sep = "__"))

#' Classify edges of two state networks as shared / A-only / B-only (+ sign switch)
#' @param edges_a Edge data.frame (node_a, node_b, r, q, sign) of state A.
#' @param edges_b Edge data.frame (node_a, node_b, r, q, sign) of state B.
#' @param label_a Label of state A.
#' @param label_b Label of state B.
#' @param sign_switch_min_r Minimum absolute correlation in both networks for an edge to count as a sign switch.
compute_rewiring <- function(edges_a, edges_b, label_a = "A", label_b = "B",
                             sign_switch_min_r = 0.5) {
  ka <- .edge_key(edges_a$node_a, edges_a$node_b)
  kb <- .edge_key(edges_b$node_a, edges_b$node_b)
  keys <- union(ka, kb)
  if (!length(keys)) {
    return(data.frame(edge_id = character(), node_a = character(), node_b = character(),
                      r_a = numeric(), r_b = numeric(), state = character(),
                      sign_switch = logical(), stringsAsFactors = FALSE))
  }
  parts <- do.call(rbind, strsplit(keys, "__", fixed = TRUE))
  u <- data.frame(edge_id = keys, node_a = parts[, 1], node_b = parts[, 2],
                  r_a = edges_a$r[match(keys, ka)], q_a = edges_a$q[match(keys, ka)],
                  r_b = edges_b$r[match(keys, kb)], q_b = edges_b$q[match(keys, kb)],
                  stringsAsFactors = FALSE)
  pa <- !is.na(u$r_a); pb <- !is.na(u$r_b)
  u$state <- ifelse(pa & pb, "shared", ifelse(pa, paste0(label_a, "_only"), paste0(label_b, "_only")))
  u$sign_switch <- pa & pb & (sign(u$r_a) != sign(u$r_b)) &
    abs(ifelse(pa, u$r_a, 0)) >= sign_switch_min_r & abs(ifelse(pb, u$r_b, 0)) >= sign_switch_min_r
  attr(u, "labels") <- c(label_a, label_b)
  u
}

#' Global and node-level rewiring summaries (PhenoRewire summarize_rewiring)
#' @param union Union edge table returned by \code{compute_rewiring()}.
#' @param label_a Label of state A.
#' @param label_b Label of state B.
summarize_rewiring <- function(union, label_a = "A", label_b = "B") {
  la <- paste0(label_a, "_only"); lb <- paste0(label_b, "_only")
  n_sh <- sum(union$state == "shared"); n_a <- sum(union$state == la); n_b <- sum(union$state == lb)
  n_ss <- sum(union$sign_switch)
  glob <- data.frame(n_shared_edges = n_sh, n_a_only_edges = n_a, n_b_only_edges = n_b,
                     n_sign_switch_edges = n_ss, total_edges_in_union = n_sh + n_a + n_b,
                     total_state_specific_edges = n_a + n_b,
                     rewiring_proportion = round((n_a + n_b) / max(1, n_sh + n_a + n_b), 4),
                     sign_switch_proportion = round(n_ss / max(1, n_sh), 4))
  names(glob)[2:3] <- c(paste0("n_", la, "_edges"), paste0("n_", lb, "_edges"))
  if (!nrow(union)) return(list(global = glob, nodes = data.frame()))
  long <- rbind(data.frame(f = union$node_a, s = union$state, ss = union$sign_switch),
                data.frame(f = union$node_b, s = union$state, ss = union$sign_switch))
  feats <- sort(unique(long$f))
  cnt <- function(cond) as.numeric(tapply(cond, factor(long$f, levels = feats), sum))
  nodes <- data.frame(feature_id = feats, shared_edges = cnt(long$s == "shared"),
                      a_only = cnt(long$s == la), b_only = cnt(long$s == lb),
                      sign_switch_edges = cnt(long$ss), stringsAsFactors = FALSE)
  nodes$rewiring_score <- nodes$a_only + nodes$b_only + nodes$sign_switch_edges
  nodes$total_union_degree <- nodes$shared_edges + nodes$a_only + nodes$b_only
  nodes$rewiring_fraction <- round((nodes$a_only + nodes$b_only) / pmax(1, nodes$total_union_degree), 4)
  nodes <- nodes[order(-nodes$rewiring_score, -nodes$sign_switch_edges, -nodes$a_only,
                       -nodes$b_only, nodes$feature_id), ]
  names(nodes)[3:4] <- c(paste0(la, "_edges"), paste0(lb, "_edges"))
  rownames(nodes) <- NULL
  list(global = glob, nodes = nodes)
}

#' PhenoRewire triage score for one state network
#' triage = 0.7 * mean percentile rank(degree centrality, betweenness, eigenvector,
#'          weighted degree) + 0.3 * (1 - q)
#' @param metrics_nodes Node table from \code{network_metrics()}.
#' @param q Numeric vector of BH-adjusted q-values.
#' @param topology_weight Weight of the topology score.
#' @param selection_weight Weight of the selection-support score.
triage_nodes <- function(metrics_nodes, q, topology_weight = 0.7, selection_weight = 0.3) {
  pr <- function(x) rank(x, ties.method = "average") / length(x)
  m <- metrics_nodes
  m$topology_priority_score <- (pr(m$degree_centrality) + pr(m$betweenness_centrality) +
                                  pr(m$eigenvector_centrality) + pr(m$weighted_degree_abs)) / 4
  m$selection_support_score <- 1 - pmin(pmax(q[m$feature_id], 0), 1)
  m$selection_support_score[is.na(m$selection_support_score)] <- 0
  m$triage_score <- topology_weight * m$topology_priority_score + selection_weight * m$selection_support_score
  m[order(-m$triage_score), ]
}
