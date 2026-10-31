# ---- MetaboRewire - Input & preprocessing ------------------------------------
# Readers for Metabolon-style exports and a PhenoRewire-compatible preprocessing
# step (presence filter -> optional median normalisation -> log transform ->
# removal of uninformative features).

.read_bom_csv <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  fileEncoding = "UTF-8-BOM", na.strings = c("", "NA", "N/A"))
}

#' Read a Metabolon dataset (one matrix / tissue)
#'
#' @param data_dir folder that holds the Metabolon csv exports.
#' @param suffix file suffix: "" for plasma, "H" for heart in this project.
#' @param tissue label stored in the returned object.
#' @param group_map named character vector mapping GROUP_NAME -> state label.
#' @return list(X, detected, meta, anno, tissue); X is features x samples
#'   (Metabolon batch-normalised, minimum-imputed values), detected is a logical
#'   matrix (TRUE where the non-imputed batch-normalised value is observed).
read_metabolon_dataset <- function(data_dir, suffix = "", tissue = "plasma",
                                   group_map = c("day 0" = "Healthy",
                                                 "day 40" = "EAM",
                                                 "day 40 Metformin" = "EAM_Met")) {
  f <- function(stem) file.path(data_dir, paste0(stem, suffix, ".csv"))
  imp  <- .read_bom_csv(f("batchnormimputed"))
  raw  <- .read_bom_csv(f("batchnorm"))
  meta <- .read_bom_csv(f("metadata"))
  anno <- .read_bom_csv(f("chemicalannotation"))

  stopifnot(identical(imp$PARENT_SAMPLE_NAME, raw$PARENT_SAMPLE_NAME))
  samples <- imp$PARENT_SAMPLE_NAME
  feat_cols <- setdiff(colnames(imp), c("PARENT_SAMPLE_NAME", "GROUP_NAME"))
  X <- t(as.matrix(imp[, feat_cols]))
  storage.mode(X) <- "double"
  colnames(X) <- samples
  D <- t(!is.na(as.matrix(raw[match(samples, raw$PARENT_SAMPLE_NAME), feat_cols])))
  colnames(D) <- samples
  rownames(D) <- rownames(X) <- feat_cols

  meta <- meta[match(samples, meta$PARENT_SAMPLE_NAME), ]
  sample_num <- as.integer(sub(".*-", "", samples))
  state <- unname(group_map[meta$GROUP_NAME])
  if (anyNA(state)) stop("Unmapped GROUP_NAME values: ",
                         paste(unique(meta$GROUP_NAME[is.na(state)]), collapse = ", "))
  meta_out <- data.frame(
    sample_id = samples,
    state = factor(state, levels = unique(unname(group_map))),
    mouse_id = as.character(meta$CLIENT_SAMPLE_ID),
    # Metabolon sample numbering arrives in two contiguous runs of 9 samples
    # (3 per group); kept as a technical 'block' for sensitivity analyses.
    block = ifelse(rank(sample_num) <= length(sample_num) / 2, "B1", "B2"),
    treatment_field = meta$TREATMENT,
    stringsAsFactors = FALSE)

  anno$CHEM_ID <- as.character(anno$CHEM_ID)
  anno <- anno[match(feat_cols, anno$CHEM_ID), ]
  anno_out <- data.frame(
    feature_id = feat_cols,
    name = ifelse(is.na(anno$CHEMICAL_NAME), feat_cols, anno$CHEMICAL_NAME),
    super_pathway = ifelse(is.na(anno$SUPER_PATHWAY) | anno$SUPER_PATHWAY == "",
                           "Unannotated", anno$SUPER_PATHWAY),
    sub_pathway = ifelse(is.na(anno$SUB_PATHWAY) | anno$SUB_PATHWAY == "",
                         "Unannotated", anno$SUB_PATHWAY),
    type = anno$TYPE, kegg = anno$KEGG, hmdb = anno$HMDB, platform = anno$PLATFORM,
    stringsAsFactors = FALSE)
  rownames(anno_out) <- feat_cols
  list(X = X, detected = D, meta = meta_out, anno = anno_out, tissue = tissue)
}

#' Column-wise median normalisation (PhenoRewire utils.median_normalize)
#' @param X Numeric matrix, features x samples.
median_normalize <- function(X) {
  med <- apply(X, 2, stats::median, na.rm = TRUE)
  target <- stats::median(med, na.rm = TRUE)
  sc <- target / med
  sc[!is.finite(sc)] <- 1
  sweep(X, 2, sc, "*")
}

#' Preprocess a features x samples matrix
#'
#' @param X numeric matrix, features x samples.
#' @param detected optional logical matrix of the same shape; when given, the
#'   presence filter uses it (Metabolon imputed data have no zeros).
#' @param min_presence_global minimum fraction of samples in which a feature is
#'   detected (PhenoRewire MIN_PRESENCE_GLOBAL).
#' @param min_nonzero_intensity values below this are set to 0 (PhenoRewire).
#' @param norm "none" or "median".
#' @param transform "log2" (for strictly positive, already scaled data such as
#'   Metabolon), "log2p1" (PhenoRewire default) or "none".
#' @param groups optional sample grouping for a group-wise presence rule.
#' @return list(X = processed matrix, kept = feature ids, log = list of counts)
preprocess_matrix <- function(X, detected = NULL, min_presence_global = 0.5,
                              min_nonzero_intensity = 0, norm = c("none", "median"),
                              transform = c("log2", "log2p1", "none"), groups = NULL) {
  norm <- match.arg(norm); transform <- match.arg(transform)
  X <- as.matrix(X); X[!is.finite(X)] <- 0
  if (min_nonzero_intensity > 0) X[X < min_nonzero_intensity] <- 0
  present <- if (is.null(detected)) X > 0 else detected[rownames(X), colnames(X)]
  # with `groups`, a feature is kept when it reaches min_presence_global in at
  # least one group (keeps presence/absence features such as drugs)
  keep <- if (is.null(groups)) rowMeans(present) >= min_presence_global else
    apply(sapply(unique(groups), function(gr) rowMeans(present[, groups == gr, drop = FALSE])), 1, max) >= min_presence_global
  n0 <- nrow(X)
  X <- X[keep, , drop = FALSE]
  raw_const <- apply(X, 1, stats::sd) <= 0
  if (norm == "median") X <- median_normalize(X)
  X <- switch(transform, log2 = log2(X), log2p1 = log2(X + 1), none = X)
  post_const <- apply(X, 1, stats::sd) <= 0 | !is.finite(apply(X, 1, stats::sd))
  drop <- raw_const | post_const
  X <- X[!drop, , drop = FALSE]
  list(X = X, kept = rownames(X),
       log = list(n_input = n0, n_presence_kept = sum(keep),
                  n_constant_removed = sum(drop), n_final = nrow(X)))
}

#' Features eligible for a state-specific network: detected in at least
#' `min_state_presence` of that state's samples. Avoids correlations driven by
#' tied imputed values (a MetaboRewire addition; set to 0 to disable).
#' @param detected Logical matrix (features x samples), TRUE where the value was measured rather than imputed.
#' @param samples Sample ids belonging to one state.
#' @param min_state_presence Minimum fraction of the state samples in which a feature must be detected.
state_presence_ok <- function(detected, samples, min_state_presence) {
  if (is.null(detected) || min_state_presence <= 0) return(NULL)
  rowMeans(detected[, samples, drop = FALSE]) >= min_state_presence
}
