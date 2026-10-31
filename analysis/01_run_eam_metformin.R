#!/usr/bin/env Rscript
# =============================================================================
# MetaboRewire - EAM / metformin study (plasma + heart, Metabolon HD4)
#
# 1. QC + preprocessing (+ technical-block adjustment)
# 2. Global multivariate structure (PCA, PERMANOVA, PERMDISP)
# 3. Differential abundance (limma, three contrasts)
# 4. PhenoRewire-faithful pairwise rewiring (3 comparisons x 2 tissues)
#    + MetaboRewire exact permutation inference and jackknife edge support
#    + sensitivity: limma-based selection, no block adjustment
# 5. Whole-metabolome threshold-free differential co-variation + pathway cross-talk
# 6. Three-state trajectory: hard-network rescue (NRI) and threshold-free rescue
# 7. Heart-plasma (inter-organ) coupling
#
# Usage:  Rscript 01_run_eam_metformin.R [config.yaml]
# =============================================================================
suppressMessages({
  library(MetaboRewire); library(limma); library(vegan); library(yaml); library(jsonlite)
})
# Arguments come from the command line (Rscript) or, in RStudio, from MR_ARGS set by RUN_ALL.R
args <- if (exists("MR_ARGS")) MR_ARGS else commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args)) args[1] else system.file("config", "eam_metformin.yaml", package = "MetaboRewire")
cfg <- yaml::read_yaml(cfg_path)
if (length(args) >= 2) cfg$data_dir <- args[2]
if (length(args) >= 3) cfg$out_dir <- args[3]
out <- cfg$out_dir
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
wcsv <- function(x, ...) utils::write.csv(x, file.path(out, ...), row.names = FALSE)
mkd <- function(...) { d <- file.path(out, ...); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
states <- unlist(cfg$states)
gmap <- unlist(cfg$group_map)
pr <- cfg$phenorewire
pr$covariate <- cfg$preprocess$covariate
key <- list()        # key numbers for the manuscript

# ---------------------------------------------------------------------------
# 1. Load, QC, preprocess
# ---------------------------------------------------------------------------
D <- list()
for (tis in names(cfg$tissues)) {
  ds <- read_metabolon_dataset(cfg$data_dir, cfg$tissues[[tis]], tis, gmap)
  st <- as.character(ds$meta$state)
  pp <- preprocess_matrix(ds$X, ds$detected, cfg$preprocess$min_presence_group,
                          norm = cfg$preprocess$norm, transform = cfg$preprocess$transform, groups = st)
  X_raw <- pp$X
  X <- if (isTRUE(cfg$preprocess$adjust_covariate) && !is.null(cfg$preprocess$covariate))
    adjust_covariate(X_raw, ds$meta[[cfg$preprocess$covariate]], st) else X_raw
  # network-eligible features: detected in >= min_state_presence of every state
  pres_state <- sapply(states, function(s) rowMeans(ds$detected[rownames(X), st == s, drop = FALSE]))
  net_ok <- apply(pres_state, 1, min) >= pr$min_state_presence
  D[[tis]] <- list(ds = ds, X = X, X_raw = X_raw, state = st, net_ok = net_ok, pp_log = pp$log,
                   classes = stats::setNames(ds$anno[rownames(X), cfg$whole_metabolome$class_column], rownames(X)))
  # QC: drug detection and metadata consistency
  met_id <- ds$anno$feature_id[tolower(ds$anno$name) == "metformin"]
  qc <- data.frame(sample_id = ds$meta$sample_id, mouse_id = ds$meta$mouse_id, state = st,
                   block = ds$meta$block, TREATMENT_field = ds$meta$treatment_field,
                   metformin_detected = if (length(met_id)) ds$detected[met_id, ds$meta$sample_id] else NA)
  qc$flag <- ifelse(qc$state == "EAM_Met" & qc$TREATMENT_field != "metformin",
                    "TREATMENT field says no treatment but sample is in metformin group", "")
  wcsv(qc, paste0("qc_samples_", tis, ".csv"))
  key[[tis]]$preprocess <- c(pp$log, n_network_eligible = sum(net_ok))
  key[[tis]]$metformin_detected <- tapply(qc$metformin_detected, qc$state, sum)
  key[[tis]]$metadata_flags <- qc$sample_id[qc$flag != ""]
  log_msg(tis, ": ", nrow(X), " features; ", sum(net_ok), " network-eligible")
}

# ---------------------------------------------------------------------------
# 2. Global multivariate structure
# ---------------------------------------------------------------------------
mv <- list()
for (tis in names(D)) {
  d <- D[[tis]]; meta <- d$ds$meta
  for (adj in c("raw", "adjusted")) {
    X <- if (adj == "raw") d$X_raw else d$X
    S <- scale(t(X))
    pc <- prcomp(S, center = FALSE)
    mv[[tis]][[adj]] <- list(scores = data.frame(meta, PC1 = pc$x[, 1], PC2 = pc$x[, 2]),
                             var = summary(pc)$importance[2, 1:5])
  }
  S <- scale(t(d$X_raw)); dst <- dist(S)
  set.seed(42)
  a <- vegan::adonis2(dst ~ state + block, data = meta, by = "margin", permutations = 9999)
  tab <- data.frame(term = rownames(a), as.data.frame(a), row.names = NULL)
  pw <- do.call(rbind, lapply(cfg$comparisons, function(cm) {
    idx <- meta$state %in% cm
    set.seed(42)
    r <- vegan::adonis2(dist(S[idx, ]) ~ state + block, data = droplevels(meta[idx, ]), by = "margin",
                        permutations = 9999)
    bd <- vegan::betadisper(dist(S[idx, ]), droplevels(meta$state[idx]))
    set.seed(42); pd <- vegan::permutest(bd, permutations = 9999)
    data.frame(tissue = tis, comparison = paste(cm, collapse = " vs "), R2_state = r$R2[1],
               F_state = r$F[1], p_state = r$`Pr(>F)`[1], R2_block = r$R2[2], p_block = r$`Pr(>F)`[2],
               permdisp_F = pd$tab$F[1], permdisp_p = pd$tab$`Pr(>F)`[1])
  }))
  pw$p_state_BH <- p.adjust(pw$p_state, "BH")
  wcsv(tab, paste0("permanova_global_", tis, ".csv")); wcsv(pw, paste0("permanova_pairwise_", tis, ".csv"))
  key[[tis]]$permanova <- tab; key[[tis]]$permanova_pairwise <- pw
}
log_msg("multivariate done")

# ---------------------------------------------------------------------------
# 3. Differential abundance (limma; state + block)
# ---------------------------------------------------------------------------
DA <- list()
for (tis in names(D)) {
  d <- D[[tis]]; meta <- d$ds$meta
  st <- factor(d$state, levels = states)
  des <- model.matrix(~ 0 + st + factor(meta$block)); colnames(des) <- c("H", "E", "M", "blockB2")
  fit <- eBayes(contrasts.fit(lmFit(d$X_raw, des),
                              makeContrasts(EAM_vs_Healthy = E - H, Met_vs_EAM = M - E,
                                            Met_vs_Healthy = M - H, levels = des)))
  res <- do.call(rbind, lapply(colnames(fit$coefficients), function(k) {
    tt <- topTable(fit, coef = k, number = Inf, sort.by = "none")
    data.frame(tissue = tis, contrast = k, feature_id = rownames(tt), log2FC = tt$logFC,
               t = tt$t, p = tt$P.Value, q = tt$adj.P.Val, stringsAsFactors = FALSE)
  }))
  res <- merge(res, d$ds$anno[, c("feature_id", "name", "super_pathway", "sub_pathway", "kegg", "hmdb")],
               by = "feature_id")
  res <- res[order(res$contrast, res$p), ]
  DA[[tis]] <- res
  wcsv(res, paste0("differential_abundance_", tis, ".csv"))
  key[[tis]]$DA <- aggregate(cbind(q05 = q < 0.05, q10 = q < 0.10, p05 = p < 0.05) ~ contrast, res, sum)
  # metabolite-level abundance rescue for disease-changed metabolites
  w <- reshape(res[, c("feature_id", "contrast", "log2FC", "q")], idvar = "feature_id",
               timevar = "contrast", direction = "wide")
  w$abundance_rescue <- 1 - w$log2FC.Met_vs_Healthy / w$log2FC.EAM_vs_Healthy
  w <- merge(w, d$ds$anno[, c("feature_id", "name", "super_pathway", "sub_pathway")], by = "feature_id")
  wcsv(w[w$q.EAM_vs_Healthy < 0.10, ], paste0("abundance_rescue_disease_metabolites_", tis, ".csv"))
}
log_msg("DA done")

# ---------------------------------------------------------------------------
# 4. PhenoRewire-faithful pairwise rewiring (+ MetaboRewire inference)
# ---------------------------------------------------------------------------
PR <- list(); pr_tab <- list()
variants <- list(primary = list(method = "rank_perm", adjusted = TRUE),
                 limma_selection = list(method = "limma", adjusted = TRUE),
                 no_block_adjustment = list(method = "rank_perm", adjusted = FALSE))
for (tis in names(D)) for (vn in names(variants)) for (cm in cfg$comparisons) {
  d <- D[[tis]]; v <- variants[[vn]]
  cfg_run <- pr; cfg_run$selection_method <- v$method
  if (!v$adjusted) cfg_run$covariate <- NULL
  od <- mkd("phenorewire", vn, tis, paste(cm, collapse = "_vs_"))
  r <- run_phenorewire(if (v$adjusted) d$X else d$X_raw, d$ds$meta, cm[1], cm[2], d$ds$anno,
                       d$ds$detected, cfg_run, od)
  nm <- paste(tis, vn, cm[1], cm[2], sep = "|")
  if (vn == "primary") PR[[nm]] <- r
  g <- if (r$status == "ok") stats::setNames(r$perm$global$observed, r$perm$global$stat) else NULL
  gp <- if (r$status == "ok") stats::setNames(r$perm$global$p_value, r$perm$global$stat) else NULL
  gm <- if (r$status == "ok") stats::setNames(r$perm$global$null_mean, r$perm$global$stat) else NULL
  pr_tab[[nm]] <- data.frame(
    tissue = tis, variant = vn, reference = cm[1], case = cm[2],
    n_tested = nrow(r$selection), applied_fdr = r$selection_meta$applied_fdr_alpha,
    n_selected = r$n_selected, n_network_nodes = r$n_network_features, status = r$status,
    edges_ref = if (!is.null(g)) g["n_a"] else NA, edges_case = if (!is.null(g)) g["n_b"] else NA,
    shared = if (!is.null(g)) g["n_shared"] else NA,
    state_specific = if (!is.null(g)) g["n_specific"] else NA,
    state_specific_null_mean = if (!is.null(g)) gm["n_specific"] else NA,
    p_state_specific = if (!is.null(g)) gp["n_specific"] else NA,
    rewiring_proportion = if (!is.null(g)) g["rewiring_proportion"] else NA,
    rewiring_proportion_null_mean = if (!is.null(g)) gm["rewiring_proportion"] else NA,
    modularity_ref = if (!is.null(g)) r$metrics_ref$summary$modularity else NA,
    modularity_case = if (!is.null(g)) r$metrics_case$summary$modularity else NA,
    DC = if (!is.null(g)) g["DC_global"] else NA, DC_null_mean = if (!is.null(g)) gm["DC_global"] else NA,
    p_DC = if (!is.null(g)) gp["DC_global"] else NA,
    coherence_ref = if (!is.null(g)) g["coherence_a"] else NA,
    coherence_case = if (!is.null(g)) g["coherence_b"] else NA,
    p_delta_coherence = if (!is.null(g)) gp["delta_coherence"] else NA,
    n_nodes_pDC_05 = if (!is.null(g)) sum(r$perm$nodes$p_DC <= 0.05) else NA,
    mean_jackknife_support = if (!is.null(g)) mean(c(r$net_ref$edges$jackknife_support,
                                                     r$net_case$edges$jackknife_support)) else NA,
    row.names = NULL)
  log_msg("phenorewire ", nm, " -> ", r$status)
}
pr_tab <- do.call(rbind, pr_tab)
wcsv(pr_tab, "phenorewire_pairwise_summary.csv")

# ---------------------------------------------------------------------------
# 5. Whole-metabolome differential co-variation + pathway cross-talk
# ---------------------------------------------------------------------------
WM <- list()
for (tis in names(D)) {
  d <- D[[tis]]; Xn <- d$X[d$net_ok, ]
  for (cm in cfg$comparisons) {
    nm <- paste(tis, cm[1], cm[2], sep = "|")
    r <- global_dc_test(Xn, d$state, cm, d$classes[rownames(Xn)], n_perm = cfg$whole_metabolome$n_perm)
    r$nodes <- merge(r$nodes, d$ds$anno[, c("feature_id", "name", "super_pathway", "sub_pathway")], by = "feature_id")
    r$nodes <- r$nodes[order(r$nodes$p_DC), ]
    od <- mkd("whole_metabolome_dc", tis)
    tag <- paste(cm, collapse = "_vs_")
    utils::write.csv(r$nodes, file.path(od, paste0("node_DC_", tag, ".csv")), row.names = FALSE)
    utils::write.csv(r$pairs, file.path(od, paste0("pathway_crosstalk_", tag, ".csv")), row.names = FALSE)
    utils::write.csv(r$global, file.path(od, paste0("global_", tag, ".csv")), row.names = FALSE)
    WM[[nm]] <- r
    log_msg("whole-metabolome DC ", nm)
  }
}
wm_tab <- do.call(rbind, lapply(names(WM), function(nm) {
  r <- WM[[nm]]; s <- strsplit(nm, "\\|")[[1]]
  data.frame(tissue = s[1], reference = s[2], case = s[3], n_features = nrow(r$nodes),
             global_DC = r$global$observed[1], global_DC_null = r$global$null_mean[1], p_global_DC = r$global$p[1],
             delta_coherence = r$global$observed[2], p_delta_coherence = r$global$p[2],
             nodes_p01 = sum(r$nodes$p_DC <= 0.01), nodes_expected_p01 = 0.01 * nrow(r$nodes),
             nodes_q10 = sum(r$nodes$q_DC <= 0.10), pathway_pairs_p05 = sum(r$pairs$p <= 0.05),
             pathway_pairs_q10 = sum(r$pairs$q <= 0.10), n_pathway_pairs = nrow(r$pairs))
}))
wcsv(wm_tab, "whole_metabolome_dc_summary.csv")

# ---------------------------------------------------------------------------
# 6. Three-state trajectory analysis
# ---------------------------------------------------------------------------
TS <- list(); ts_tab <- list(); tsd_tab <- list()
for (tis in names(D)) {
  d <- D[[tis]]; Xn <- d$X[d$net_ok, ]
  # (a) hard networks on Kruskal-Wallis-selected features
  sel <- select_features(Xn, factor(d$state, levels = states), "multigroup",
                         n_perm = cfg$three_state$n_perm_selection, fdr_alpha = pr$selection_fdr,
                         adaptive = TRUE, min_selected = pr$selection_min_features,
                         alpha_ceiling = pr$selection_fdr_ceiling, min_features_hard_stop = 0)
  fs <- sel$feature_id[sel$selected]
  ks_fdr <- attr(sel, "meta")$applied_fdr_alpha
  p <- do.call(network_params, pr$network)
  od <- mkd("three_state", tis)
  sel <- merge(sel, d$ds$anno[, c("feature_id", "name", "super_pathway", "sub_pathway")], by = "feature_id")
  utils::write.csv(sel[order(sel$q), ], file.path(od, "kruskal_selection.csv"), row.names = FALSE)
  hard <- NULL
  if (length(fs) >= 5) {
    tr <- three_state_rewiring(Xn[fs, ], d$state, states, p)
    rp <- rescue_permutation_test(Xn[fs, ], d$state, states, tr, p, n_perm = cfg$three_state$n_perm)
    tr$edges$name_a <- d$ds$anno[tr$edges$node_a, "name"]; tr$edges$name_b <- d$ds$anno[tr$edges$node_b, "name"]
    tr$edges$class_a <- d$classes[tr$edges$node_a]; tr$edges$class_b <- d$classes[tr$edges$node_b]
    utils::write.csv(tr$edges, file.path(od, "three_state_edges.csv"), row.names = FALSE)
    nodes <- merge(tr$nodes, d$ds$anno[, c("feature_id", "name", "super_pathway", "sub_pathway")], by = "feature_id")
    utils::write.csv(nodes[order(-nodes$disease_perturbed), ], file.path(od, "three_state_nodes.csv"), row.names = FALSE)
    utils::write.csv(rp$table, file.path(od, "rescue_permutation.csv"), row.names = FALSE)
    ct <- if (nrow(tr$edges)) pathway_crosstalk(tr$edges, d$classes[fs]) else data.frame()
    utils::write.csv(ct, file.path(od, "three_state_pathway_crosstalk.csv"), row.names = FALSE)
    hard <- list(tr = tr, rp = rp, n_features = length(fs), applied_fdr = ks_fdr)
    ts_tab[[tis]] <- data.frame(tissue = tis, n_features = length(fs),
                                applied_fdr = ks_fdr,
                                abs_r = tr$thresholds$abs_r, fdr = tr$thresholds$fdr,
                                edges_Healthy = tr$n_edges_state[1], edges_EAM = tr$n_edges_state[2],
                                edges_EAM_Met = tr$n_edges_state[3], t(tr$counts),
                                NRI_null_mean = rp$table$null_mean[rp$table$stat == "NRI"],
                                p_NRI = rp$table$p_upper[rp$table$stat == "NRI"],
                                p_treatment_specific = rp$table$p_upper[rp$table$stat == "treatment_specific"],
                                row.names = NULL)
  }
  # (b) threshold-free rescue (whole network-eligible metabolome)
  soft <- three_state_dc(Xn, d$state, states, d$classes[rownames(Xn)], n_perm = cfg$three_state$n_perm)
  soft$nodes <- merge(soft$nodes, d$ds$anno[, c("feature_id", "name", "super_pathway", "sub_pathway")], by = "feature_id")
  utils::write.csv(soft$global, file.path(od, "threshold_free_rescue_global.csv"), row.names = FALSE)
  utils::write.csv(soft$pairs[order(soft$pairs$p_reversal), ], file.path(od, "pathway_crosstalk_rescue.csv"), row.names = FALSE)
  utils::write.csv(soft$nodes, file.path(od, "node_coherence_by_state.csv"), row.names = FALSE)
  tsd_tab[[tis]] <- data.frame(tissue = tis, soft$global, row.names = NULL)
  TS[[tis]] <- list(hard = hard, soft = soft)
  log_msg("three-state ", tis)
}
wcsv(do.call(rbind, ts_tab), "three_state_hard_summary.csv")
wcsv(do.call(rbind, tsd_tab), "three_state_threshold_free_summary.csv")

# ---------------------------------------------------------------------------
# 7. Heart-plasma (inter-organ) coupling, paired by animal
# ---------------------------------------------------------------------------
CP <- NULL
if (all(c("plasma", "heart") %in% names(D))) {
  dh <- D$heart; dp <- D$plasma
  to_mouse <- function(dd) { X <- dd$X; colnames(X) <- dd$ds$meta$mouse_id; X }
  Xh <- to_mouse(dh)[dh$net_ok, ]; Xp <- to_mouse(dp)[dp$net_ok, ]
  st_m <- stats::setNames(dh$state, dh$ds$meta$mouse_id)
  stopifnot(all(st_m[dp$ds$meta$mouse_id] == dp$state))
  cc <- cross_compartment_coupling(Xh, Xp, st_m, n_perm = cfg$coupling$n_perm)
  cc <- merge(cc, dh$ds$anno[, c("feature_id", "name", "super_pathway", "sub_pathway")], by = "feature_id")
  od <- mkd("heart_plasma_coupling")
  utils::write.csv(cc[order(cc$p_pooled), ], file.path(od, "same_metabolite_coupling.csv"), row.names = FALSE)
  # state-wise mean coupling and paired permutation test of its change
  ids <- intersect(colnames(Xh), colnames(Xp)); f <- cc$feature_id
  A <- Xh[f, ids]; B <- Xp[f, ids]; s <- st_m[ids]
  mean_coupling <- function(lab) sapply(states, function(k) {
    rr <- rowSums(rank_rows(A[, lab == k]) * rank_rows(B[, lab == k])); mean(rr) })
  # standardise within state so that permutation removes only the dependence structure
  As <- standardize_within(A, s); Bs <- standardize_within(B, s)
  mc_obs <- sapply(states, function(k) mean(rowSums(rank_rows(As[, s == k]) * rank_rows(Bs[, s == k]))))
  set.seed(42)
  null <- t(replicate(cfg$coupling$n_perm, { ps <- sample(s)
    sapply(states, function(k) mean(rowSums(rank_rows(As[, ps == k]) * rank_rows(Bs[, ps == k])))) }))
  dif <- function(m) { m <- unname(m); c(EAM_minus_Healthy = m[2] - m[1], Met_minus_EAM = m[3] - m[2], Met_minus_Healthy = m[3] - m[1]) }
  d_obs <- dif(mc_obs); d_null <- t(apply(null, 1, dif))
  coupling_global <- data.frame(contrast = names(d_obs), observed = d_obs,
                                p = sapply(seq_along(d_obs), function(j) mean(abs(d_null[, j]) >= abs(d_obs[j]) - 1e-12)))
  bycls <- do.call(rbind, lapply(split(cc, cc$super_pathway), function(z)
    data.frame(super_pathway = z$super_pathway[1], n = nrow(z), t(colMeans(z[, paste0("rho_", states)])),
               mean_within_pooled = mean(z$rho_within_pooled))))
  utils::write.csv(coupling_global, file.path(od, "coupling_change_permutation.csv"), row.names = FALSE)
  utils::write.csv(data.frame(state = states, mean_coupling = mc_obs), file.path(od, "mean_coupling_by_state.csv"), row.names = FALSE)
  utils::write.csv(bycls, file.path(od, "coupling_by_pathway.csv"), row.names = FALSE)
  # bipartite inter-organ rewiring on disease-associated metabolites (limma q < 0.10 in each tissue)
  bip <- list()
  for (cm in cfg$comparisons) {
    ctr <- switch(paste(cm, collapse = "_"), Healthy_EAM = "EAM_vs_Healthy", EAM_EAM_Met = "Met_vs_EAM",
                  Healthy_EAM_Met = "Met_vs_Healthy")
    fh <- head(DA$heart$feature_id[DA$heart$contrast == ctr & DA$heart$q < 0.10 & DA$heart$feature_id %in% rownames(Xh)], 60)
    fp <- head(DA$plasma$feature_id[DA$plasma$contrast == ctr & DA$plasma$q < 0.10 & DA$plasma$feature_id %in% rownames(Xp)], 60)
    if (length(fh) < 3 || length(fp) < 3) next
    m <- s %in% cm
    br <- bipartite_rewiring(Xh[fh, ids[m]], Xp[fp, ids[m]], s[m], list(abs_r = 0.8, fdr = 0.2))
    if (!is.null(br$edges)) {
      br$edges$name_heart <- dh$ds$anno[br$edges$feature_A, "name"]
      br$edges$name_plasma <- dp$ds$anno[br$edges$feature_B, "name"]
      utils::write.csv(br$edges, file.path(od, paste0("bipartite_edges_", paste(cm, collapse = "_vs_"), ".csv")), row.names = FALSE)
    }
    bip[[paste(cm, collapse = "_vs_")]] <- data.frame(comparison = paste(cm, collapse = " vs "),
      n_heart = length(fh), n_plasma = length(fp), br$global, row.names = NULL)
  }
  bip_tab <- do.call(rbind, bip)
  if (!is.null(bip_tab)) utils::write.csv(bip_tab, file.path(od, "bipartite_rewiring_summary.csv"), row.names = FALSE)
  CP <- list(cc = cc, coupling_global = coupling_global, mean_coupling = mc_obs, bycls = bycls, bip = bip_tab)
  log_msg("coupling done")
}

# ---------------------------------------------------------------------------
# Save everything for figures / manuscript
# ---------------------------------------------------------------------------
saveRDS(list(D = lapply(D, function(d) d[c("X", "X_raw", "state", "net_ok", "classes", "pp_log")]),
             anno = lapply(D, function(d) d$ds$anno), meta = lapply(D, function(d) d$ds$meta),
             mv = mv, DA = DA, PR = lapply(PR, function(r) r[setdiff(names(r), c("X"))]),
             pr_tab = pr_tab, WM = WM, wm_tab = wm_tab, TS = TS, CP = CP, key = key, cfg = cfg),
        file.path(out, "metaborewire_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
log_msg("all done -> ", out)
