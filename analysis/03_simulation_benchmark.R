#!/usr/bin/env Rscript
# =============================================================================
# MetaboRewire - Simulation benchmark: calibration and power of rewiring
# statistics at small sample sizes.
#
# Scenarios (p = 60 metabolites, two states, n per state in {4, 6, 10, 20}):
#   null    : identical dependence structure in both states (3 modules of 15
#             metabolites, within-module rho = 0.7, 15 independent metabolites);
#             20 metabolites additionally shifted by 1 SD in state B (abundance
#             change without rewiring).
#   rewired : as null, but in state B module 1 decouples (rho = 0) and a new
#             module (rho = 0.7) forms among 10 previously independent
#             metabolites -> 25 truly rewired metabolites.
#   block   : null structure plus a balanced technical block (half of the
#             samples of each state) shifting all metabolites by 1 SD
#             (tests spurious edges and the effect of covariate adjustment).
#
# Methods: PhenoRewire hard-threshold statistics (rewiring proportion, number of
# state-specific edges) reported naively and calibrated by permutation, and the
# MetaboRewire threshold-free differential connectivity (DC) test.
# Usage: Rscript 03_simulation_benchmark.R [n_rep] [n_perm] [out_dir] [cores]
# =============================================================================
suppressMessages(library(MetaboRewire))
# Arguments come from the command line (Rscript) or, in RStudio, from MR_ARGS set by RUN_ALL.R
args <- if (exists("MR_ARGS")) MR_ARGS else commandArgs(trailingOnly = TRUE)
n_rep  <- if (length(args) >= 1) as.integer(args[1]) else 100
n_perm <- if (length(args) >= 2) as.integer(args[2]) else 200
out_dir <- if (length(args) >= 3) args[3] else "results/simulation"
cores <- if (length(args) >= 4) as.integer(args[4]) else max(1, parallel::detectCores() - 1)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

P <- 60
mod <- list(1:15, 16:30, 31:45)
indep <- 46:60
new_mod <- 46:55
true_rewired <- c(mod[[1]], new_mod)

make_sigma <- function(modules, rho = 0.7) {
  S <- diag(P)
  for (m in modules) { S[m, m] <- rho; diag(S)[m] <- 1 }
  S
}
S_A <- make_sigma(mod)
S_B_alt <- make_sigma(list(mod[[2]], mod[[3]], new_mod))
chol_A <- chol(S_A); chol_B <- chol(S_B_alt)
rmv <- function(n, ch) matrix(stats::rnorm(n * P), n) %*% ch

simulate <- function(scenario, n) {
  A <- rmv(n, chol_A)
  B <- if (scenario == "rewired") rmv(n, chol_B) else rmv(n, chol_A)
  B[, 1:20] <- B[, 1:20] + 1                      # abundance shift, no rewiring
  X <- t(rbind(A, B)); rownames(X) <- sprintf("m%02d", 1:P)
  g <- rep(c("A", "B"), each = n)
  blk <- rep(rep(c("b1", "b2"), each = n / 2), 2)
  if (scenario == "block") X[, blk == "b2"] <- X[, blk == "b2"] + 1
  list(X = X, g = g, blk = blk)
}

auc <- function(score, truth) {
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

true_edges <- function(modules) {
  k <- c(); for (m in modules) { cm <- utils::combn(sort(m), 2); k <- c(k, paste(sprintf("m%02d", cm[1, ]), sprintf("m%02d", cm[2, ]))) }
  k
}
TE_A <- true_edges(mod)

one_rep <- function(scenario, n, rep) {
  set.seed(1e5 * match(scenario, c("null", "rewired", "block")) + 1e3 * n + rep)
  d <- simulate(scenario, n)
  p <- network_params()
  out <- list()
  for (adj in if (scenario == "block") c(FALSE, TRUE) else FALSE) {
    X <- if (adj) adjust_covariate(d$X, d$blk, d$g) else d$X
    pt <- rewiring_permutation_test(X, d$g, p, n_perm = n_perm, seed = rep, levels = c("A", "B"))
    gl <- stats::setNames(pt$global$observed, pt$global$stat)
    pv <- stats::setNames(pt$global$p_value, pt$global$stat)
    truth <- rownames(X) %in% sprintf("m%02d", true_rewired)
    nd <- pt$nodes[match(rownames(X), pt$nodes$feature_id), ]
    # false-edge rate in the state-A network (spurious edges from the block effect)
    na <- correlation_network(X[, d$g == "A"], p)
    fe <- if (nrow(na$edges)) mean(!(paste(pmin(na$edges$node_a, na$edges$node_b),
                                           pmax(na$edges$node_a, na$edges$node_b)) %in% TE_A)) else NA
    out[[length(out) + 1]] <- data.frame(
      scenario = scenario, n = n, rep = rep, block_adjusted = adj,
      n_edges_A = gl["n_a"], n_edges_B = gl["n_b"], n_specific = gl["n_specific"],
      rewiring_proportion = gl["rewiring_proportion"],
      naive_call = gl["rewiring_proportion"] >= 0.5,
      p_hard = pv["n_specific"], p_jaccard = pv["jaccard_distance"], p_DC = pv["DC_global"],
      auc_hard = if (scenario == "rewired") auc(nd$rewiring_score, truth) else NA,
      auc_DC = if (scenario == "rewired") auc(nd$DC, truth) else NA,
      false_edge_rate_A = fe, row.names = NULL)
  }
  do.call(rbind, out)
}

grid <- expand.grid(rep = seq_len(n_rep), n = c(4, 6, 10, 20),
                    scenario = c("null", "rewired", "block"), stringsAsFactors = FALSE)
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(grid)), function(i)
  tryCatch(one_rep(grid$scenario[i], grid$n[i], grid$rep[i]), error = function(e) NULL),
  mc.cores = cores)
res <- do.call(rbind, res)
message("simulation time: ", format(Sys.time() - t0))
utils::write.csv(res, file.path(out_dir, "simulation_replicates.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(res, list(res$scenario, res$n, res$block_adjusted), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], n = d$n[1], block_adjusted = d$block_adjusted[1], reps = nrow(d),
             mean_rewiring_proportion = mean(d$rewiring_proportion),
             naive_call_rate = mean(d$naive_call),
             reject_hard_perm = mean(d$p_hard <= 0.05),
             reject_jaccard_perm = mean(d$p_jaccard <= 0.05),
             reject_DC_perm = mean(d$p_DC <= 0.05),
             mean_auc_hard = mean(d$auc_hard), mean_auc_DC = mean(d$auc_DC),
             mean_false_edge_rate_A = mean(d$false_edge_rate_A, na.rm = TRUE))))
summ <- summ[order(summ$scenario, summ$block_adjusted, summ$n), ]
utils::write.csv(summ, file.path(out_dir, "simulation_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
