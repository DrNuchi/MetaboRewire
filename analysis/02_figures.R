#!/usr/bin/env Rscript
# =============================================================================
# MetaboRewire - publication figures (PDF + 300-dpi PNG)
# Usage: Rscript 02_figures.R <results_dir> <simulation_dir> <fig_dir>
# =============================================================================
suppressMessages({ library(MetaboRewire); library(ggplot2); library(patchwork); library(ggrepel) })
# Arguments come from the command line (Rscript) or, in RStudio, from MR_ARGS set by RUN_ALL.R
args <- if (exists("MR_ARGS")) MR_ARGS else commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else "results/eam_metformin"
sim_dir <- if (length(args) >= 2) args[2] else "results/simulation"
fig_dir <- if (length(args) >= 3) args[3] else "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
R <- readRDS(file.path(res_dir, "metaborewire_results.rds"))

# ---- design tokens (validated categorical palette; ink tokens for text) -----
col_state <- c(Healthy = "#2a78d6", EAM = "#eb6834", EAM_Met = "#1baf7a")
lab_state <- c(Healthy = "Healthy (d0)", EAM = "EAM (d40)", EAM_Met = "EAM + metformin (d40)")
shp_state <- c(Healthy = 21, EAM = 24, EAM_Met = 22)
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"; neutral <- "#f0efec"
div_lo <- "#2a78d6"; div_hi <- "#e34948"
theme_mr <- theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2),
        panel.grid.major = element_line(colour = grid_col, linewidth = 0.25),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 8.5),
        plot.tag = element_text(face = "bold", size = 10), legend.key.size = unit(3, "mm"),
        strip.text = element_text(face = "bold", colour = ink), plot.background = element_rect(fill = "white", colour = NA))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf")
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}
tis_lab <- c(plasma = "Plasma", heart = "Heart")
cmp_lab <- function(a, b) paste0(sub("EAM_Met", "Met", b), " vs ", sub("EAM_Met", "Met", a))

# =============================================================================
# Figure 2 - study structure: PCA, DA counts, PERMANOVA
# =============================================================================
pca_panel <- function(tis) {
  s <- R$mv[[tis]]$adjusted$scores; v <- R$mv[[tis]]$adjusted$var
  ggplot(s, aes(PC1, PC2, fill = state, shape = state)) +
    stat_ellipse(aes(colour = state), type = "t", level = 0.68, linewidth = 0.3, show.legend = FALSE) +
    geom_point(size = 2.2, colour = "white", stroke = 0.4) +
    scale_fill_manual(values = col_state, labels = lab_state, name = NULL) +
    scale_colour_manual(values = col_state) +
    scale_shape_manual(values = shp_state, labels = lab_state, name = NULL) +
    labs(title = paste0(tis_lab[tis], " (block-adjusted)"),
         x = sprintf("PC1 (%.1f%%)", 100 * v[1]), y = sprintf("PC2 (%.1f%%)", 100 * v[2])) + theme_mr
}
da <- do.call(rbind, lapply(names(R$DA), function(t) {
  d <- R$DA[[t]]
  do.call(rbind, lapply(split(d, d$contrast), function(z) data.frame(
    tissue = tis_lab[t], contrast = z$contrast[1],
    direction = c("up", "down"), n = c(sum(z$q < 0.10 & z$log2FC > 0), sum(z$q < 0.10 & z$log2FC < 0))))) }))
da$contrast <- factor(da$contrast, levels = c("EAM_vs_Healthy", "Met_vs_EAM", "Met_vs_Healthy"),
                      labels = c("EAM vs Healthy", "Met vs EAM", "Met vs Healthy"))
da$tissue <- factor(da$tissue, levels = c("Plasma", "Heart"))
da$signed <- ifelse(da$direction == "up", da$n, -da$n)
p_da <- ggplot(da, aes(contrast, signed, fill = direction)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n, vjust = ifelse(signed >= 0, -0.3, 1.2)), size = 2.4, colour = ink2) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  facet_wrap(~tissue) + scale_fill_manual(values = c(up = div_hi, down = div_lo),
                                          labels = c(up = "Higher", down = "Lower"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(title = "Differential abundance (limma, q < 0.10)", x = NULL, y = "Metabolites") +
  theme_mr + theme(axis.text.x = element_text(angle = 25, hjust = 1))
pm <- do.call(rbind, lapply(names(R$key), function(t) { z <- R$key[[t]]$permanova_pairwise; z$tissue <- tis_lab[t]; z }))
pm$comparison <- factor(gsub("EAM_Met", "Met", pm$comparison), levels = c("Healthy vs EAM", "EAM vs Met", "Healthy vs Met"))
pm$tissue <- factor(pm$tissue, levels = c("Plasma", "Heart"))
p_pm <- ggplot(pm, aes(comparison, R2_state, fill = tissue)) +
  geom_col(position = position_dodge(0.7), width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(p_state_BH < 0.001, "<.001", sub("^0", "", sprintf("%.3f", p_state_BH)))),
            position = position_dodge(0.7), vjust = -0.4, size = 1.9, colour = ink2) +
  scale_fill_manual(values = c(Plasma = "#6250d6", Heart = "#e87ba4"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "PERMANOVA, state | block (BH q)", x = NULL, y = expression(R^2~"(state)")) + theme_mr
fig2 <- (pca_panel("plasma") + pca_panel("heart") + plot_layout(guides = "collect")) /
  (p_da + p_pm + plot_layout(widths = c(1.4, 1))) + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
save_fig(fig2, "Fig3_study_overview", 175, 140)

# =============================================================================
# Figure 3 - calibration: simulation + observed vs permutation null
# =============================================================================
sim_file <- file.path(sim_dir, "simulation_summary.csv")
if (file.exists(sim_file)) {
  ss <- read.csv(sim_file); sr <- read.csv(file.path(sim_dir, "simulation_replicates.csv"))
  sr2 <- sr[sr$scenario %in% c("null", "rewired") & !sr$block_adjusted, ]
  sr2$scenario <- factor(sr2$scenario, levels = c("null", "rewired"), labels = c("Null (no rewiring)", "True rewiring"))
  pa <- ggplot(sr2, aes(factor(n), rewiring_proportion, fill = scenario)) +
    geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.3, position = position_dodge(0.7), colour = ink2) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = ink2, linewidth = 0.3) +
    scale_fill_manual(values = c("#c3c2b7", "#eb6834"), name = NULL) +
    labs(title = "Naive rewiring proportion", x = "Samples per state", y = "State-specific / union edges") + theme_mr
  s2 <- ss[ss$scenario %in% c("null", "rewired") & !ss$block_adjusted, ]
  lp <- rbind(data.frame(n = s2$n, scenario = s2$scenario, method = "Naive (proportion >= 0.5)", rate = s2$naive_call_rate),
              data.frame(n = s2$n, scenario = s2$scenario, method = "Hard edges + permutation", rate = s2$reject_hard_perm),
              data.frame(n = s2$n, scenario = s2$scenario, method = "Threshold-free DC + permutation", rate = s2$reject_DC_perm))
  lp$scenario <- factor(lp$scenario, levels = c("null", "rewired"), labels = c("Type I error (null)", "Power (rewired)"))
  col_m <- c("Naive (proportion >= 0.5)" = "#52514e", "Hard edges + permutation" = "#eb6834",
             "Threshold-free DC + permutation" = "#2a78d6")
  pb <- ggplot(lp, aes(n, rate, colour = method, shape = method)) +
    geom_hline(data = data.frame(scenario = factor("Type I error (null)", levels = levels(lp$scenario)), y = 0.05),
               aes(yintercept = y), linetype = 2, colour = ink2, linewidth = 0.3) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) + facet_wrap(~scenario) +
    scale_x_continuous(breaks = c(4, 6, 10, 20), trans = "log2") + scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = col_m, name = NULL) + scale_shape_manual(values = c("Naive (proportion >= 0.5)" = 17, "Hard edges + permutation" = 16, "Threshold-free DC + permutation" = 15), name = NULL) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    labs(title = "Global test calibration and power", x = "Samples per state", y = "Rejection rate") + theme_mr
  rw <- sr[sr$scenario == "rewired", ]
  au <- rbind(data.frame(n = sort(unique(rw$n)), method = "Hard rewiring score", auc = tapply(rw$auc_hard, rw$n, mean, na.rm = TRUE)),
              data.frame(n = sort(unique(rw$n)), method = "Threshold-free DC", auc = tapply(rw$auc_DC, rw$n, mean, na.rm = TRUE)))
  pc <- ggplot(au, aes(n, auc, colour = method, shape = method)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = ink2, linewidth = 0.3) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    scale_x_continuous(breaks = c(4, 6, 10, 20), trans = "log2") + scale_y_continuous(limits = c(0.3, 1)) +
    scale_colour_manual(values = c("Hard rewiring score" = "#eb6834", "Threshold-free DC" = "#2a78d6"), name = NULL) +
    scale_shape_manual(values = c("Hard rewiring score" = 16, "Threshold-free DC" = 15), name = NULL) +
    labs(title = "Node ranking (AUROC)", x = "Samples per state", y = "AUROC, rewired metabolites") + theme_mr
  sb <- ss[ss$scenario == "block", ]
  sb$cond <- ifelse(sb$block_adjusted, "Block, adjusted", "Block, unadjusted")
  sn <- ss[ss$scenario == "null", ]; sn$cond <- "No block (reference)"
  sb <- rbind(sb, sn)
  pd <- ggplot(sb, aes(n, mean_false_edge_rate_A, colour = cond, shape = cond)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    scale_x_continuous(breaks = c(4, 6, 10, 20), trans = "log2") + scale_y_log10() +
    scale_colour_manual(values = c("Block, unadjusted" = "#e34948", "Block, adjusted" = "#1baf7a", "No block (reference)" = "#52514e"), name = NULL) +
    scale_shape_manual(values = c("Block, unadjusted" = 16, "Block, adjusted" = 15, "No block (reference)" = 4), name = NULL) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    labs(title = "Technical block: spurious edges", x = "Samples per state", y = "False-edge fraction (log)") + theme_mr
} else { pa <- pb <- pc <- pd <- ggplot() + theme_void() }
pt <- R$pr_tab[R$pr_tab$variant == "primary" & R$pr_tab$status == "ok", ]
pt$label <- paste0(tis_lab[pt$tissue], ": ", cmp_lab(pt$reference, pt$case))
pt$label <- factor(pt$label, levels = rev(pt$label))
pe_has <- pt$edges_ref + pt$edges_case > 0
obs_null <- rbind(data.frame(label = pt$label[pe_has], stat = "Rewiring proportion", observed = pt$rewiring_proportion[pe_has],
                             null = pt$rewiring_proportion_null_mean[pe_has], p = pt$p_state_specific[pe_has]),
                  data.frame(label = pt$label, stat = "Differential connectivity (DC)", observed = pt$DC,
                             null = pt$DC_null_mean, p = pt$p_DC))
obs_null$stat <- factor(obs_null$stat, levels = c("Rewiring proportion", "Differential connectivity (DC)"))
pe <- ggplot(obs_null, aes(y = label)) +
  geom_segment(aes(x = null, xend = observed, yend = label), colour = ink2, linewidth = 0.4) +
  geom_point(aes(x = null), shape = 21, fill = "white", colour = ink2, size = 2) +
  geom_point(aes(x = observed), shape = 21, fill = "#eb6834", colour = "white", size = 2.3) +
  geom_text(aes(x = pmax(observed, null), label = sprintf("p=%.3f", p)), hjust = -0.25, size = 2.1, colour = ink2) +
  facet_wrap(~stat, scales = "free_x") + scale_x_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  labs(title = "EAM data: observed (filled) vs permutation null mean (open); p for state-specific edges / DC", x = NULL, y = NULL) + theme_mr
fig3 <- (pa + pb + plot_layout(widths = c(1, 1.6))) / (pc + pd) / pe + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
save_fig(fig3, "Fig2_calibration", 175, 215)

# =============================================================================
# Figure 4 - pathway cross-talk: change in mean |rho| between super-pathways
# =============================================================================
abbr <- c("Amino Acid" = "AA", "Carbohydrate" = "Carb", "Cofactors and Vitamins" = "Cof/Vit", "Energy" = "Energy",
          "Lipid" = "Lipid", "Nucleotide" = "Nucl", "Partially Characterized Molecules" = "PCM", "Peptide" = "Pept",
          "Xenobiotics" = "Xeno", "Unannotated" = "Unann")
ct <- do.call(rbind, lapply(names(R$WM), function(nm) {
  s <- strsplit(nm, "\\|")[[1]]; p <- R$WM[[nm]]$pairs
  p2 <- p; p2$class_a <- p$class_b; p2$class_b <- p$class_a
  z <- unique(rbind(p, p2)); z$tissue <- tis_lab[s[1]]; z$cmp <- cmp_lab(s[2], s[3]); z }))
ct <- ct[ct$cmp != "Met vs Healthy", ]
ct$cmp <- factor(ct$cmp, levels = c("EAM vs Healthy", "Met vs EAM"))
ct$tissue <- factor(ct$tissue, levels = c("Plasma", "Heart"))
ct$a <- factor(abbr[ct$class_a], levels = abbr); ct$b <- factor(abbr[ct$class_b], levels = abbr)
ct <- ct[!is.na(ct$a) & !is.na(ct$b) & as.integer(ct$a) >= as.integer(ct$b), ]
lim <- max(abs(ct$delta_coherence))
fig4 <- ggplot(ct, aes(a, b, fill = delta_coherence)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))), size = 3, colour = ink, vjust = 0.75) +
  facet_grid(tissue ~ cmp) +
  scale_fill_gradient2(low = div_lo, mid = neutral, high = div_hi, midpoint = 0, limits = c(-lim, lim),
                       name = expression(Delta~"mean |"*rho*"|")) +
  labs(x = NULL, y = NULL, title = "Pathway cross-talk rewiring (* p < 0.05, ** p < 0.01; exact permutation)") +
  coord_equal() + theme_mr + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank())
save_fig(fig4, "Fig4_pathway_crosstalk", 175, 150)

# =============================================================================
# Figure 5 - three-state rescue: network vs abundance, edge classes
# =============================================================================
rs <- do.call(rbind, lapply(names(R$TS), function(t) {
  g <- R$TS[[t]]$soft$global; nl <- R$TS[[t]]$soft$null
  rbind(data.frame(tissue = tis_lab[t], level = "Network coherence", null = nl[, 1], obs = g$observed_cor[1], p = g$p_rescue[1]),
        data.frame(tissue = tis_lab[t], level = "Abundance", null = nl[, 2], obs = g$observed_cor[2], p = g$p_rescue[2])) }))
rs$tissue <- factor(rs$tissue, levels = c("Plasma", "Heart"))
rs$level <- factor(rs$level, levels = c("Network coherence", "Abundance"))
lab5 <- unique(rs[, c("tissue", "level", "obs", "p")])
p5a <- ggplot(rs, aes(null)) + geom_histogram(bins = 30, fill = "#c3c2b7", colour = "white", linewidth = 0.2) +
  geom_vline(data = lab5, aes(xintercept = obs), colour = "#eb6834", linewidth = 0.7) +
  geom_text(data = lab5, aes(x = Inf, y = Inf, label = sprintf("observed (line)\np = %.3f", p)), hjust = 1.05, vjust = 1.3,
            size = 2.2, colour = ink, lineheight = 0.9) +
  facet_grid(tissue ~ level, scales = "free") +
  labs(title = "Rescue: cor(EAM change, metformin change) across metabolites vs exchangeable null",
       x = "Correlation across metabolites (more negative = more reversal)", y = "Permutations") + theme_mr
ec <- do.call(rbind, lapply(names(R$TS), function(t) { h <- R$TS[[t]]$hard; if (is.null(h)) return(NULL)
  tb <- h$rp$table; cc <- h$tr$counts
  cls <- c("disease_gained_rescued", "disease_gained_persistent", "disease_lost_restored", "disease_lost_persistent",
           "treatment_emergent", "treatment_disrupted", "conserved", "sign_switch")
  nm <- colMeans(h$rp$null[, cls, drop = FALSE], na.rm = TRUE)
  data.frame(tissue = tis_lab[t], class = cls, observed = as.numeric(cc[cls]), null = as.numeric(nm)) }))
if (!is.null(ec)) {
  ec$class <- factor(ec$class, levels = rev(unique(ec$class)),
                     labels = rev(c("Gained in EAM, reversed", "Gained in EAM, persistent", "Lost in EAM, restored",
                                    "Lost in EAM, persistent", "Treatment-emergent", "Treatment-disrupted",
                                    "Conserved", "Sign switch")))
  ec$tissue <- factor(ec$tissue, levels = c("Plasma", "Heart"))
  p5b <- ggplot(ec, aes(y = class)) +
    geom_col(aes(x = observed), fill = "#2a78d6", width = 0.6) +
    geom_point(aes(x = null), shape = 21, fill = "white", colour = ink, size = 1.8) +
    facet_wrap(~tissue, scales = "free_x") +
    labs(title = "Edge classes vs null mean (o)", x = "Edges", y = NULL) + theme_mr
} else p5b <- ggplot() + theme_void()
pp <- do.call(rbind, lapply(names(R$TS), function(t) { z <- R$TS[[t]]$soft$pairs; z$tissue <- tis_lab[t]; z }))
pp$tissue <- factor(pp$tissue, levels = c("Plasma", "Heart"))
pp$lab <- paste(abbr[pp$class_a], abbr[pp$class_b], sep = "-")
p5c <- ggplot(pp, aes(delta_disease, delta_treatment)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) + geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_abline(slope = -1, intercept = 0, linetype = 2, colour = ink2, linewidth = 0.3) +
  geom_point(aes(fill = p_reversal < 0.05), shape = 21, colour = "white", size = 2) +
  geom_text_repel(data = do.call(rbind, lapply(split(pp, pp$tissue), function(z) head(z[order(z$p_reversal), ][z$p_reversal[order(z$p_reversal)] < 0.05, ], 4))), aes(label = lab), size = 2, colour = ink, max.overlaps = 20,
                  segment.size = 0.2) +
  scale_fill_manual(values = c("FALSE" = "#c3c2b7", "TRUE" = "#eb6834"), labels = c("p >= 0.05", "p < 0.05"),
                    name = "Reversal") +
  facet_wrap(~tissue) +
  labs(title = "Pathway-pair coherence change",
       x = expression(Delta~"mean |"*rho*"| (EAM - Healthy)"), y = expression(Delta~"mean |"*rho*"| (Met - EAM)")) + theme_mr
fig5 <- p5a / (p5b + p5c + plot_layout(widths = c(1, 1.3))) + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
save_fig(fig5, "Fig5_three_state_rescue", 175, 170)

# =============================================================================
# Figure 6 - heart-plasma coupling
# =============================================================================
if (!is.null(R$CP)) {
  cc <- R$CP$cc
  long <- do.call(rbind, lapply(names(col_state), function(s)
    data.frame(state = s, rho = cc[[paste0("rho_", s)]], super = cc$super_pathway)))
  long$state <- factor(long$state, levels = names(col_state))
  cg <- R$CP$coupling_global
  p6a <- ggplot(long, aes(state, rho, fill = state)) +
    geom_violin(colour = NA, alpha = 0.35, scale = "width", width = 0.8) +
    geom_boxplot(width = 0.18, outlier.shape = NA, colour = ink2, linewidth = 0.3) +
    scale_fill_manual(values = col_state, guide = "none") + scale_x_discrete(labels = lab_state) +
    labs(title = "Heart-plasma coupling (same metabolite)",
         subtitle = paste(sprintf("%s: %+.3f (p = %.2f)", gsub("_", " ", cg$contrast), cg$observed, cg$p), collapse = "\n"),
         x = NULL, y = expression("Spearman"~rho~"(heart vs plasma, n = 6 mice)")) +
    theme_mr + theme(plot.subtitle = element_text(size = 6.5, colour = ink2), axis.text.x = element_text(angle = 15, hjust = 1))
  top <- head(cc[order(cc$p_pooled), ], 12)
  top$name <- factor(top$name, levels = rev(top$name))
  tl <- do.call(rbind, lapply(names(col_state), function(s) data.frame(name = top$name, state = s, rho = top[[paste0("rho_", s)]])))
  tl$state <- factor(tl$state, levels = names(col_state))
  p6b <- ggplot(tl, aes(rho, name, fill = state, shape = state, group = state)) +
    geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3) +
    geom_point(size = 2, colour = "white", stroke = 0.3, position = position_dodge(width = 0.65)) +
    scale_fill_manual(values = col_state, labels = lab_state, name = NULL) +
    scale_shape_manual(values = shp_state, labels = lab_state, name = NULL) +
    labs(title = "Top coupled metabolites (all q = 0.009)", x = expression(rho), y = NULL) + theme_mr
  fig6 <- (p6a + p6b + plot_layout(widths = c(1, 1.2))) + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
  save_fig(fig6, "Fig6_heart_plasma_coupling", 175, 85)
}

# =============================================================================
# Supplementary - PhenoRewire rewiring networks (primary pairwise runs)
# =============================================================================
net_plot <- function(r, title) {
  u <- r$union; nodes <- r$nodes
  if (!nrow(u)) return(ggplot() + annotate("text", 0, 0, label = "no edges", colour = ink2) + theme_void() + ggtitle(title))
  g <- igraph::graph_from_data_frame(u[, c("node_a", "node_b")], directed = FALSE)
  set.seed(1); L <- igraph::layout_with_fr(g)
  xy <- data.frame(id = igraph::V(g)$name, x = L[, 1], y = L[, 2])
  e <- cbind(u, xy[match(u$node_a, xy$id), c("x", "y")], setNames(xy[match(u$node_b, xy$id), c("x", "y")], c("xend", "yend")))
  la <- attr(u, "labels")
  e$cls <- ifelse(e$state == "shared", "Shared", ifelse(e$state == paste0(la[1], "_only"), paste(la[1], "only"), paste(la[2], "only")))
  xy$score <- nodes$rewiring_score[match(xy$id, nodes$feature_id)]
  xy$name <- nodes$name[match(xy$id, nodes$feature_id)]
  ecol <- stats::setNames(c(col_state[la[1]], col_state[la[2]], "#8a8a85"), c(paste(la[1], "only"), paste(la[2], "only"), "Shared"))
  ggplot() + geom_segment(data = e, aes(x, y, xend = xend, yend = yend, colour = cls), linewidth = 0.5) +
    geom_point(data = xy, aes(x, y, size = score), shape = 21, fill = "white", colour = ink2, stroke = 0.4) +
    geom_text_repel(data = if (nrow(xy) > 25) xy[xy$score >= 2, ] else xy, aes(x, y, label = substr(name, 1, 28)),
                    size = 1.8, colour = ink, max.overlaps = 30, segment.size = 0.15) +
    scale_colour_manual(values = ecol, name = NULL) + scale_size(range = c(1, 3.5), guide = "none") +
    ggtitle(title) + theme_void(base_size = 7) + theme(plot.title = element_text(face = "bold", size = 8),
                                                        legend.position = "bottom")
}
nets <- R$PR[sapply(R$PR, function(r) r$status == "ok")]
nets <- nets[sapply(nets, function(r) nrow(r$union) > 0)]
pl <- lapply(names(nets), function(nm) { s <- strsplit(nm, "\\|")[[1]]
  net_plot(nets[[nm]], paste0(tis_lab[s[1]], ": ", cmp_lab(s[3], s[4]))) })
if (length(pl)) save_fig(wrap_plots(pl, ncol = 2) + plot_annotation(tag_levels = "a"), "FigS1_phenorewire_networks", 175, 75 * ceiling(length(pl) / 2))
message("figures written to ", fig_dir)
