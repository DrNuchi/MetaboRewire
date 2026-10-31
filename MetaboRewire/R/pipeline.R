# ---- MetaboRewire - Export and one-call pairwise pipeline --------------------

#' Write a rewiring network (union graph) to GraphML for Cytoscape
#' @param union Union edge table returned by \code{compute_rewiring()}.
#' @param node_table Node attribute table with a \code{feature_id} column.
#' @param path Output GraphML file path.
#' @param label_a Label of state A.
#' @param label_b Label of state B.
export_graphml <- function(union, node_table, path, label_a, label_b) {
  if (!nrow(union)) return(invisible(NULL))
  e <- data.frame(from = union$node_a, to = union$node_b, state = union$state,
                  r_ref = ifelse(is.na(union$r_a), 0, union$r_a),
                  r_case = ifelse(is.na(union$r_b), 0, union$r_b),
                  sign_switch = union$sign_switch, stringsAsFactors = FALSE)
  v <- node_table[node_table$feature_id %in% unique(c(e$from, e$to)), , drop = FALSE]
  v <- v[!duplicated(v$feature_id), ]
  names(v)[names(v) == "feature_id"] <- "name"
  for (cn in names(v)) if (is.logical(v[[cn]])) v[[cn]] <- as.integer(v[[cn]])
  for (cn in names(v)) if (is.numeric(v[[cn]])) v[[cn]][!is.finite(v[[cn]])] <- 0
  for (cn in names(v)) if (is.character(v[[cn]])) v[[cn]][is.na(v[[cn]])] <- ""
  g <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = v)
  igraph::write_graph(g, path, format = "graphml")
  invisible(path)
}

.wcsv <- function(x, path) utils::write.csv(x, path, row.names = FALSE)

#' Run one PhenoRewire-style phenotype comparison (+ MetaboRewire inference)
#'
#' @param X features x samples (preprocessed, all samples)
#' @param meta data.frame with sample_id and state
#' @param ref,case state labels
#' @param anno feature annotation (feature_id, name, super_pathway, sub_pathway)
#' @param detected optional detection matrix for within-state presence filtering
#' @param cfg list of parameters (see inst/config/eam_metformin.yaml)
#' @param outdir output folder
run_phenorewire <- function(X, meta, ref, case, anno, detected = NULL, cfg, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  keep_s <- meta$state %in% c(ref, case)
  m <- meta[keep_s, ]; Xs <- X[, m$sample_id, drop = FALSE]
  y <- as.integer(m$state == case)
  # comparison-level presence filter (PhenoRewire MIN_PRESENCE_GLOBAL per run)
  if (!is.null(detected)) {
    pres <- rowMeans(detected[rownames(Xs), m$sample_id, drop = FALSE])
    Xs <- Xs[pres >= cfg$min_presence_run, , drop = FALSE]
  }
  Xs <- Xs[apply(Xs, 1, stats::sd) > 0, , drop = FALSE]

  # 1. feature selection: PhenoRewire rank-permutation (default) or limma
  if (identical(cfg$selection_method, "limma")) {
    cov <- if (!is.null(cfg$covariate)) m[[cfg$covariate]] else NULL
    lm_res <- limma_select(Xs, y, cov)
    q <- stats::p.adjust(lm_res$p, "BH")
    asel <- adaptive_fdr_selection(q, cfg$selection_fdr, cfg$adaptive_selection,
                                   cfg$selection_min_features, cfg$selection_fdr_ceiling)
    sel <- data.frame(feature_id = rownames(Xs), stat = lm_res$stat, p_empirical = lm_res$p,
                      q = q, selected = asel$selected,
                      mean_ref = rowMeans(Xs[, y == 0, drop = FALSE]),
                      mean_case = rowMeans(Xs[, y == 1, drop = FALSE]), stringsAsFactors = FALSE)
    sel$log2FC_case_vs_ref <- sel$mean_case - sel$mean_ref
    attr(sel, "meta") <- c(asel$meta, list(n_features_total = nrow(Xs), mode = "limma"))
  } else {
    sel <- select_features(Xs, y, "phenotype", n_perm = cfg$n_perm_selection,
                           fdr_alpha = cfg$selection_fdr, adaptive = cfg$adaptive_selection,
                           min_selected = cfg$selection_min_features,
                           alpha_ceiling = cfg$selection_fdr_ceiling, seed = cfg$seed)
  }
  selmeta <- attr(sel, "meta")
  sel <- merge(sel, anno[, c("feature_id", "name", "super_pathway", "sub_pathway")],
               by = "feature_id", all.x = TRUE)
  sel <- sel[order(sel$q, sel$p_empirical), ]
  .wcsv(sel, file.path(outdir, "phenotype_features_all.csv"))
  fs <- sel$feature_id[sel$selected]
  # MetaboRewire: network features must be detected in >= min_state_presence of BOTH states
  if (!is.null(detected) && cfg$min_state_presence > 0) {
    ok <- rowMeans(detected[fs, m$sample_id[y == 0], drop = FALSE]) >= cfg$min_state_presence &
          rowMeans(detected[fs, m$sample_id[y == 1], drop = FALSE]) >= cfg$min_state_presence
    fs_net <- fs[ok]
  } else fs_net <- fs
  sel_sel <- sel[sel$feature_id %in% fs, ]
  sel_sel$in_network_set <- sel_sel$feature_id %in% fs_net
  .wcsv(sel_sel, file.path(outdir, "phenotype_features_selected.csv"))
  res <- list(selection = sel, selection_meta = selmeta, n_selected = length(fs),
              n_network_features = length(fs_net))
  if (length(fs_net) < 3) {
    res$status <- "too_few_features"
    jsonlite::write_json(list(status = res$status, n_selected = length(fs)),
                         file.path(outdir, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
    return(res)
  }
  p <- do.call(network_params, cfg$network)

  # 2. state networks
  Xr <- Xs[fs_net, y == 0, drop = FALSE]; Xc <- Xs[fs_net, y == 1, drop = FALSE]
  na <- correlation_network(Xr, p); nb <- correlation_network(Xc, p)
  qv <- stats::setNames(sel$q, sel$feature_id)
  ma <- network_metrics(fs_net, na$edges, cfg$seed); mb <- network_metrics(fs_net, nb$edges, cfg$seed)
  na$edges$jackknife_support <- jackknife_edge_support(Xr, na$edges, na$thresholds$abs_r)
  nb$edges$jackknife_support <- jackknife_edge_support(Xc, nb$edges, nb$thresholds$abs_r)
  for (s in list(list(ref, na, ma), list(case, nb, mb))) {
    d <- file.path(outdir, paste0("network_", s[[1]])); dir.create(d, showWarnings = FALSE)
    .wcsv(s[[2]]$edges, file.path(d, "edges.csv"))
    .wcsv(triage_nodes(s[[3]]$nodes, qv), file.path(d, "priority_features.csv"))
    .wcsv(cbind(s[[3]]$summary, as.data.frame(s[[2]]$thresholds[c("abs_r", "fdr")])),
          file.path(d, "network_summary.csv"))
  }

  # 3. rewiring
  u <- compute_rewiring(na$edges, nb$edges, ref, case, p$sign_switch_min_r)
  rs <- summarize_rewiring(u, ref, case)
  # 4. permutation inference (MetaboRewire)
  pt <- rewiring_permutation_test(Xs[fs_net, ], ifelse(y == 1, case, ref), p,
                                  n_perm = cfg$n_perm_rewiring, seed = cfg$seed,
                                  levels = c(ref, case))
  nodes <- merge(if (nrow(rs$nodes)) rs$nodes else data.frame(feature_id = character()),
                 pt$nodes[, c("feature_id", "null_mean", "p_perm", "q_perm", "DC", "DC_null_mean", "p_DC", "q_DC")],
                 by = "feature_id", all = TRUE)
  for (cn in c("shared_edges", "sign_switch_edges", "rewiring_score"))
    nodes[[cn]] <- if (cn %in% names(nodes)) ifelse(is.na(nodes[[cn]]), 0, nodes[[cn]]) else 0
  nodes <- merge(nodes, sel[, c("feature_id", "name", "super_pathway", "sub_pathway", "q",
                                "log2FC_case_vs_ref")], by = "feature_id", all.x = TRUE)
  nodes <- nodes[order(-nodes$rewiring_score, nodes$p_DC), ]
  rw <- file.path(outdir, "rewiring"); dir.create(rw, showWarnings = FALSE)
  .wcsv(u, file.path(rw, "edges_union_with_stats.csv"))
  .wcsv(rs$global, file.path(rw, "rewiring_summary.csv"))
  .wcsv(nodes, file.path(rw, "priority_rewired_nodes.csv"))
  .wcsv(pt$global, file.path(rw, "rewiring_permutation_global.csv"))
  vt <- merge(data.frame(feature_id = fs_net), nodes, all.x = TRUE)
  vt$community_ref <- ma$nodes$community_id[match(vt$feature_id, ma$nodes$feature_id)]
  vt$community_case <- mb$nodes$community_id[match(vt$feature_id, mb$nodes$feature_id)]
  export_graphml(u, vt, file.path(rw, "rewiring_network_phenotype.graphml"), ref, case)

  summ <- list(comparison = paste(ref, "vs", case), n_samples_ref = sum(y == 0),
               n_samples_case = sum(y == 1), n_features_tested = nrow(Xs),
               selection = selmeta, n_selected = length(fs), n_network_features = length(fs_net),
               network_ref = c(as.list(ma$summary), na$thresholds),
               network_case = c(as.list(mb$summary), nb$thresholds),
               rewiring = as.list(rs$global),
               permutation = list(exact = pt$exact, n_perm = pt$n_perm,
                                  global = pt$global))
  jsonlite::write_json(summ, file.path(outdir, "summary.json"), auto_unbox = TRUE,
                       pretty = TRUE, digits = 6)
  c(res, list(status = "ok", features = fs_net, net_ref = na, net_case = nb, metrics_ref = ma,
              metrics_case = mb, union = u, rewiring = rs, perm = pt, nodes = nodes,
              X = Xs[fs_net, ], y = y, summary = summ))
}
