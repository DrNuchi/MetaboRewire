test_that("rank-based association equals Spearman correlation", {
  set.seed(1)
  X <- matrix(rnorm(80), 10); y <- rnorm(8)
  expect_equal(spearman_assoc(X, y), apply(X, 1, cor, y = y, method = "spearman"))
})

test_that("exact permutation enumerates all labelings and gives exact p", {
  X <- rbind(c(1, 2, 3, 10, 11, 12), c(3, 1, 2, 2, 3, 1))
  r <- perm_test_binary(X, c(0, 0, 0, 1, 1, 1), n_perm = "exact")
  expect_true(r$exact); expect_equal(r$n_perm, choose(6, 3))
  expect_equal(r$p[1], 2 / 20)          # observed labeling and its mirror
  expect_equal(r$r[1], cor(X[1, ], c(0, 0, 0, 1, 1, 1), method = "spearman"))
})

test_that("adaptive FDR relaxation follows the PhenoRewire rule", {
  q <- c(0.03, 0.08, 0.09, 0.15, 0.15, 0.15, 0.5)
  a <- adaptive_fdr_selection(q, 0.05, TRUE, 5, 0.2)
  expect_equal(a$meta$applied_fdr_alpha, 0.2); expect_equal(sum(a$selected), 6)
  b <- adaptive_fdr_selection(q, 0.05, FALSE, 5, 0.2)
  expect_equal(sum(b$selected), 1)
})

test_that("threshold grid matches PhenoRewire defaults", {
  G <- threshold_grid(0.7, 0.05, 0.5, 0.2)
  expect_equal(sort(unique(G$abs_r), decreasing = TRUE), c(0.8, 0.75, 0.7, 0.65, 0.6, 0.5))
  expect_equal(sort(unique(G$fdr)), c(0.01, 0.05, 0.1, 0.2))
  expect_equal(G$abs_r[1], 0.8); expect_equal(G$fdr[1], 0.01)   # same iteration order as Python
})

test_that("Spearman p-values use the t approximation and r = 1 gives p = 0", {
  X <- rbind(1:6, c(1, 3, 2, 5, 4, 6), 6:1)
  cp <- spearman_cor_p(X)
  r <- cp$R[1, 2]; t <- r * sqrt(4 / (1 - r^2))
  expect_equal(cp$P[1, 2], 2 * pt(-abs(t), 4))
  expect_equal(cp$P[1, 3], 0)
})

test_that("edge classification: shared, state-specific and sign switch", {
  ea <- data.frame(node_a = c("a", "a", "b"), node_b = c("b", "c", "c"), r = c(0.9, 0.8, -0.7), q = 0.01)
  eb <- data.frame(node_a = c("b", "c", "c"), node_b = c("a", "d", "b"), r = c(0.95, 0.9, 0.75), q = 0.01)
  u <- compute_rewiring(ea, eb, "A", "B", 0.5)
  expect_equal(sum(u$state == "shared"), 2)
  expect_equal(sum(u$state == "A_only"), 1); expect_equal(sum(u$state == "B_only"), 1)
  expect_equal(sum(u$sign_switch), 1)
  s <- summarize_rewiring(u, "A", "B")
  expect_equal(s$global$rewiring_proportion, 0.5)
  expect_equal(s$nodes$rewiring_score[s$nodes$feature_id == "c"], 3)
})

test_that("within-state standardisation and covariate adjustment", {
  set.seed(2)
  X <- matrix(rnorm(60), 5); g <- rep(c("A", "B"), each = 6)
  Z <- standardize_within(X, g)
  expect_equal(unname(rowMeans(Z[, g == "B"])), rep(0, 5), tolerance = 1e-12)
  blk <- rep(c("b1", "b2"), 6)
  Xb <- X; Xb[, blk == "b2"] <- Xb[, blk == "b2"] + 5
  Xa <- adjust_covariate(Xb, blk, g)
  expect_lt(max(abs(rowMeans(Xa[, blk == "b2"]) - rowMeans(Xa[, blk == "b1"]) -
                      (rowMeans(X[, blk == "b2"]) - rowMeans(X[, blk == "b1"])))), 2)
  expect_equal(adjust_covariate(X, rep("b1", 12), g), X)
})

test_that("permutation test is exact and calibrated in a small null example", {
  set.seed(3)
  X <- matrix(rnorm(15 * 8), 15); rownames(X) <- paste0("f", 1:15)
  g <- rep(c("A", "B"), each = 4)
  pt <- rewiring_permutation_test(X, g, network_params(), n_perm = "exact", levels = c("A", "B"))
  expect_true(pt$exact); expect_equal(pt$n_perm, choose(8, 4))
  expect_true(all(pt$global$p_value > 0 & pt$global$p_value <= 1))
})

test_that("three-state classes are assigned from presence patterns", {
  ks <- list(list(key = c(1, 2, 3), r = c(0.9, 0.9, 0.9)),
             list(key = c(1, 4, 5), r = c(0.9, 0.9, 0.9)),
             list(key = c(1, 2, 5, 6), r = c(0.9, 0.9, -0.9, 0.9)))
  ed <- classify_three_state_edges(ks)
  cls <- setNames(ed$class, ed$key)
  expect_equal(unname(cls[c("1", "2", "3", "4", "5", "6")]),
               c("conserved", "disease_lost_restored", "disease_lost_persistent",
                 "disease_gained_rescued", "sign_switch", "treatment_emergent"))
})
