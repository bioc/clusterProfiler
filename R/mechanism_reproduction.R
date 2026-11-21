## Reproduction rate: does a single-cell-line level-B conclusion have a second dataset behind it?
##
## Why this layer. Among the 11 datasets with known ground truth, the two rows that stopped at
## level B (`single-cell-line`) -- GR in airway and in GSE286890 -- **both point at the correct
## driver**. So B does not mean wrong, B means "there is only one piece of evidence". But "only
## one" and "several point the same way" read completely differently, and the old output used one
## sentence for both.
##
## The criterion comes from `mechanism-benchmark-spec.md` Section 4: **concordance**, not
## "agrees with the pathway direction". The latter is a sign error -- an activating driver that is
## knocked down
## should push the pathway **down**, the **opposite** of the pathway direction. That scores a true
## driver as 0 reproduction, and the truer the driver the closer to 0 (measured, HIF1A x hypoxia
## pathway: old criterion 0.00, new criterion 1.00, majority_sign = -1).
##
## Direction is decided without `pathway_evidence()`'s own es: the sign of the members' mean shift
## in the perturbation profile is used directly, otherwise the same statistic would calibrate its
## own wording (phase one Section 7.10 recorded one case of self-reference).

#' Concordance of a pathway shift across independent perturbation datasets
#'
#' For each (driver, pathway) pair, takes the mean log2 fold change of the
#' pathway members in every available perturbation dataset of that driver, and
#' reports how many of them move the pathway the same way. The majority sign is
#' the **mode of the signs**, not the sign of the sum: one dataset with a very
#' large effect must not overturn two smaller ones that agree with each other.
#' Datasets with no measured members are dropped before counting, so `k` is the
#' number that actually contributed, not the number declared in `meta`.
#'
#' This is the input [pathway_evidence()] needs to decide how strongly a
#' single-cell-line level-B row may be worded.
#'
#' @param profiles Numeric matrix, rows are genes (symbols), columns are
#'   perturbation dataset identifiers.
#' @param meta Data frame with columns `dataset` and `driver`. Optional
#'   `biosample` (for the independent-sample count) and `logFC` (the knockdown
#'   strength, used by `max_log2fc`).
#' @param pathways Named list: pathway name to a character vector of member
#'   genes. [pathway_evidence()] builds this from the enrichment result.
#' @param drivers Character vector of drivers to score. Defaults to every
#'   driver present in `meta`.
#' @param min_datasets Integer. Fewer datasets than this cannot show concordance
#'   (one dataset is trivially concordant), so the pair is reported as
#'   `not-computable`.
#' @param min_biosamples Integer. A row only counts as reproduced when the
#'   supporting datasets come from at least this many distinct known biosamples.
#'   Missing or blank biosamples do not establish independent replication.
#' @param min_genes Integer. Datasets with fewer finite pathway members than this
#'   are excluded before counting `k` or computing concordance. A pair with fewer
#'   than `min_datasets` remaining is reported as `not-computable`.
#' @param max_log2fc Optional numeric. When given and `meta` has a `logFC`
#'   column, datasets weaker than this knockdown are dropped before scoring.
#'
#' @return A data frame with one row per (driver, pathway): `driver`, `pathway`,
#'   `k` (datasets with usable measurements), `n_biosamples`, `n_genes` (members
#'   found in `profiles`), `concordance`, `majority_sign`, and `grade` which is
#'   one of `reproduces`, `consistent`, `not-reproduced` or `not-computable`.
#'
#' @examples
#' set.seed(1)
#' profiles <- matrix(rnorm(60), nrow = 20,
#'                    dimnames = list(paste0("g", 1:20), c("d1", "d2", "d3")))
#' profiles[1:5, ] <- -1
#' meta <- data.frame(dataset = c("d1", "d2", "d3"),
#'                    driver = "TF1", biosample = c("A", "B", "B"))
#' mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:5)))
#'
#' @export
mechanism_reproduction <- function(profiles, meta, pathways, drivers = NULL,
                                   min_datasets = 2L, min_biosamples = 2L,
                                   min_genes = 3L, max_log2fc = NULL) {
    if (!is.matrix(profiles) || !is.numeric(profiles))
        stop("`profiles` must be a numeric matrix.", call. = FALSE)
    if (is.null(rownames(profiles)) || is.null(colnames(profiles)))
        stop("`profiles` must have gene row names and dataset column names.",
             call. = FALSE)
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (!all(c("dataset", "driver") %in% names(meta)))
        stop("`meta` must have `dataset` and `driver` columns.", call. = FALSE)
    if (!is.list(pathways) || is.null(names(pathways)))
        stop("`pathways` must be a named list of gene vectors.", call. = FALSE)
    if (anyDuplicated(rownames(profiles)) || anyDuplicated(colnames(profiles)))
        stop("`profiles` must have unique gene and dataset names.", call. = FALSE)
    meta$dataset <- as.character(meta$dataset)
    meta$driver <- as.character(meta$driver)
    meta <- unique(meta)
    if (anyDuplicated(meta$dataset))
        stop("`meta` has inconsistent duplicate dataset identities.", call. = FALSE)
    meta <- meta[meta$dataset %in% colnames(profiles), , drop = FALSE]
    ## Knockdown-strength filter. Kept before the subsetting, otherwise `k` would count weak
    ## knockdowns -- measured, 3 of 12 drivers (ARNT / ESR1 / NR3C1) have every experiment weaker
    ## than -1, scoring 4 / 6 / 1 datasets unfiltered and zero filtered. The caller has to see that
    ## difference, so `max_log2fc` is a parameter rather than hard-coded.
    if (!is.null(max_log2fc) && "logFC" %in% names(meta)) {
        keep <- !is.na(meta$logFC) & meta$logFC <= max_log2fc
        meta <- meta[keep, , drop = FALSE]
    }
    if (is.null(drivers)) drivers <- unique(meta$driver)
    drivers <- as.character(drivers)
    has_bio <- "biosample" %in% names(meta)

    rows <- vector("list", length(pathways) * length(drivers))
    n <- 0L
    for (p in names(pathways)) {
        members <- intersect(as.character(pathways[[p]]), rownames(profiles))
        for (dr in drivers) {
            n <- n + 1L
            sel <- meta[meta$driver == dr, , drop = FALSE]
            sub <- profiles[members, sel$dataset, drop = FALSE]
            finite <- is.finite(sub)
            sub[!finite] <- NA_real_
            shift <- colMeans(sub, na.rm = TRUE)
            ## The member floor applies to each contributing dataset, not matrix rows.
            ok <- colSums(finite) >= min_genes & is.finite(shift)
            shift <- shift[ok]
            sel <- sel[ok, , drop = FALSE]
            k <- nrow(sel)
            bio <- if (has_bio) trimws(as.character(sel$biosample)) else character(0)
            nb <- if (has_bio) length(unique(bio[!is.na(bio) & nzchar(bio)])) else NA_integer_
            if (k < min_datasets) {
                rows[[n]] <- data.frame(
                    driver = dr, pathway = p, k = k, n_biosamples = nb,
                    n_genes = length(members), concordance = NA_real_,
                    majority_sign = NA_integer_, grade = "not-computable",
                    stringsAsFactors = FALSE)
                next
            }
            ## Majority = **the mode of the signs**, not the sign of sum(shift).
            ## The latter flips when effect sizes are uneven: shift = c(-0.1, -0.1, +5) is 2/3
            ## negative, yet the sum is positive, so "mostly the same way" becomes not-reproduced.
            s <- sign(shift)
            npos <- sum(s > 0); nneg <- sum(s < 0)
            maj <- if (npos == nneg) NA_integer_
                   else if (npos > nneg) 1L else -1L
            if (is.na(maj)) {
                conc <- NA_real_
                grade <- "not-computable"
            } else {
                conc <- mean(s == maj)
                grade <- if (conc >= 0.9 && k >= 3L &&
                             !is.na(nb) && nb >= min_biosamples) "reproduces"
                         else if (conc >= 0.5) "consistent"
                         else "not-reproduced"
            }
            rows[[n]] <- data.frame(
                driver = dr, pathway = p, k = k, n_biosamples = nb,
                n_genes = length(members), concordance = conc,
                majority_sign = maj, grade = grade, stringsAsFactors = FALSE)
        }
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    attr(out, "min_datasets") <- as.integer(min_datasets)
    out
}
