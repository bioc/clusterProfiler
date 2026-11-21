## Preflight gate: before any attribution runs, decide whether this data can support one.
##
## Why this layer exists. Five of eleven datasets failed, and the reasons were **visible before
## running**: two datasets held no GO term at all, one held only 23 up-regulated genes, one driver
## had a single dataset in the inventory. The old approach ran anyway and printed a row of
## `C unsupported` -- the reader could not tell "not found" from "could not be looked for".
##
## The criteria come from the applicability boundary already recorded in
## `pathway-mechanism-design.md` Section 7.9 (enough differential genes, enough pathway members, a driver
## with >=2 cell lines), **not** from back-inferring the results of this comparison.
## That also means Section 7.9 gave only qualitative statements ("23 is clearly too few, 351 clearly
## enough"), so the three thresholds below are **numbers taken between two known points, never
## measured**. The gate therefore reports a **tier plus the numbers behind it**, not a pass/fail
## verdict; a caller can move the thresholds for its own data.

## Up-regulated gene count: counted with log2FC > 1.
##
## The strict definition needs padj < 0.05 and |log2FC| >= 1, but that requires the caller to pass
## padj, and this layer only receives `signed`. Using log2FC > 1 loosens the boundary a little
## (hypoxia: 351 strict -> about 400), which does not affect the "enough or not" judgement, but the
## documentation has to say clearly that it is not the same number.
.mechanism_preflight_counts <- function(x, signed) {
    d <- if (is.data.frame(x)) x else as.data.frame(x)
    ## Prefer the genes the enrichment analysis **actually consumed** (`enrichResult`'s `@gene`).
    ## Falling back to `signed > 1` inflates the count several-fold (hypoxia: 351 strict, 881 on the
    ## fallback), so the "too few up-regulated genes" criterion would almost never fire.
    n_up <- if (inherits(x, "enrichResult") && length(methods::slot(x, "gene"))) {
        length(methods::slot(x, "gene"))
    } else if (!is.null(signed) && length(signed)) {
        sum(signed > 1, na.rm = TRUE)
    } else NA_integer_
    list(n_pathways = nrow(d), n_up = n_up, d = d)
}

#' Screen input completeness and reference coverage before driver comparison
#'
#' Run this before [pathway_evidence()]. It returns a heuristic tier, not evidence
#' of causation: `attribution` means the checked input and coverage floors pass;
#' `description` means a floor failed or a required check is unknown; `none`
#' means there are no pathways to describe.
#'
#' The boundaries are the ones already written down in the pathway-mechanism
#' design (differential genes, pathway size, driver coverage), not thresholds
#' fitted to a benchmark. They are still crude: the design notes record that 23
#' up-regulated genes are clearly too few and 351 clearly enough, so the cut
#' points below sit between two measured points rather than on a measured one.
#' The numbers behind the tier are returned with it, so a caller can move the
#' cut points and see which datasets change tier.
#'
#' @param x An enrichment result, or a data.frame of pathways.
#' @param signed Named numeric of log fold changes for every gene measured.
#'   If neither these measurements nor the enrichment input gene count are
#'   available, the gene-count check is unknown and the tier is `description`.
#' @param reference Optional `mechanismReference`; used to count candidate cell
#'   lines and background programmes. Without it, candidate coverage is queried
#'   from the inventory.
#' @param drivers Optional candidate drivers. At least one must meet the coverage
#'   floor for `attribution`. Omitted candidates or failed inventory lookups leave
#'   coverage unknown and cap the tier at `description`. A driver absent from a
#'   supplied reference has known zero coverage in that reference.
#' @param min_up Up-regulated genes required before attribution is offered.
#' @param min_members Median pathway size required. This is the pathway's own
#'   size, not how many of its members the reference happens to cover: the
#'   12-member floor recorded in the design belongs to the perturbation test's
#'   p-value, not to this descriptive input screen.
#' @param min_cell_lines Cell lines a driver needs before attribution is offered.
#' @param min_background Comparable reference programmes required for the
#'   specificity percentile to mean anything.
#' @return A `mechanismPreflight` list with `tier`, `reasons`, `metrics` and
#'   `counts`.
#' @seealso [pathway_evidence()], [mechanism_report()]
#' @export
mechanism_preflight <- function(x, signed = NULL, reference = NULL, drivers = NULL,
                                min_up = 50L, min_members = 5L,
                                min_cell_lines = 2L, min_background = 10L) {
    cnt <- .mechanism_preflight_counts(x, signed)
    reasons <- character(0)
    if (!cnt$n_pathways) {
        return(.mechanism_preflight_out(
            "none", c("no-pathways: the enrichment result holds no pathways, so there is nothing to describe or attribute."),
            cnt, NA_real_, NA_integer_, NA_integer_))
    }

    ## Pathway size is the **pathway's own** size, not its covered member count.
    ##
    ## The first version used covered members with a threshold of 8, and that blocked three
    ## successful datasets too (their median coverage was only 5-7). The 8 was copied from Section 7.9,
    ## where "12 or more members" means `test_gene_set` cannot produce a p-value -- **that is a
    ## boundary of the old statistic, not of this one**. A shape correlation is still computable on
    ## a 5-member pathway. Carrying an old boundary over to a new statistic is a category error.
    sz <- vapply(seq_len(nrow(cnt$d)), function(i) length(.pathway_members(cnt$d[i, ])), 0L)
    med_members <- stats::median(sz, na.rm = TRUE)

    n_cell_lines <- NA_integer_
    ## Initialise explicitly: `exists()` inherits from parent environments, so a same-named global
    ## would pollute the output. Control flow uses the real state (NA = not computed yet), not
    ## whether a variable exists.
    n_cell_lines_min <- NA_integer_
    per <- vapply(drivers, function(g) {
        supplied <- inherits(reference, "mechanismReference")
        inv <- if (supplied) reference$meta else
            tryCatch(.mechanism_inventory(g, NULL), error = function(e) NULL)
        if (is.null(inv) || !"perturbed_gene" %in% names(inv)) return(NA_integer_)
        d <- inv[inv$perturbed_gene %in% .mechanism_driver_names(g), , drop = FALSE]
        if (!nrow(d)) {
            if (supplied) return(0L)
            status <- attr(inv, "query_status")
            ## Only a successful empty lookup establishes zero inventory coverage.
            if (!is.null(status) && all(c("driver", "status") %in% names(status)) &&
                any(status$driver %in% g & status$status %in% "empty")) return(0L)
            return(NA_integer_)
        }
        if (!"cell_line" %in% names(d)) return(NA_integer_)
        cell <- trimws(as.character(d$cell_line))
        if (anyNA(cell) || any(!nzchar(cell))) return(NA_integer_)
        length(unique(cell))
    }, 0L)
    ## A covered candidate can pass even when other candidates cannot be assessed.
    known <- per[!is.na(per)]
    if (length(known)) {
        n_cell_lines <- max(known)
        n_cell_lines_min <- min(known)
    }
    if (!length(per) || (anyNA(per) &&
        (is.na(n_cell_lines) || n_cell_lines < min_cell_lines))) {
        reasons <- c(reasons,
            "coverage-unknown: candidate cell-line coverage could not be checked; attribution is not approved.")
    }
    if (is.na(cnt$n_up)) {
        reasons <- c(reasons,
            "up-genes-unknown: no enrichment input gene count or signed measurements were supplied; the gene-count floor is unchecked.")
    }

    if (!is.na(cnt$n_up) && cnt$n_up < min_up) {
        reasons <- c(reasons, sprintf(
            "few-up-genes: %d genes have log2FC > 1; below %d, pathways are too small to carry a driver claim.",
            cnt$n_up, min_up))
    }
    if (is.finite(med_members) && med_members < min_members) {
        reasons <- c(reasons, sprintf(
            "small-pathways: the median pathway has %d members; below %d there is nothing to test.",
            as.integer(med_members), min_members))
    }
    if (!is.na(n_cell_lines) && n_cell_lines < min_cell_lines) {
        reasons <- c(reasons, sprintf(
            "single-cell-line: the best known candidate coverage is %d cell line(s); below the heuristic floor of %d.",
            n_cell_lines, min_cell_lines))
    }
    if (!is.null(reference) && inherits(reference, "mechanismReference")) {
        n_prog <- length(unique(stats::na.omit(as.character(reference$meta$perturbed_gene))))
        if (n_prog < min_background) {
            reasons <- c(reasons, sprintf(
                "thin-reference: only %d comparable programmes; below %d the specificity percentile is unreliable.",
                n_prog, min_background))
        }
    }

    tier <- if (length(reasons)) "description" else "attribution"
    out <- .mechanism_preflight_out(tier, reasons, cnt, med_members, n_cell_lines, cnt$n_up)
    if (!is.na(n_cell_lines_min)) out$metrics$min_cell_lines <- n_cell_lines_min
    out
}

.mechanism_preflight_out <- function(tier, reasons, cnt, med_members, n_cell_lines, n_up) {
    out <- structure(list(
        tier = tier, reasons = reasons,
        metrics = data.frame(
            n_pathways = cnt$n_pathways, n_up = n_up,
            median_members = med_members, n_cell_lines = n_cell_lines,
            stringsAsFactors = FALSE),
        counts = cnt[c("n_pathways", "n_up")]),
        class = "mechanismPreflight")
    out
}

#' @rdname mechanism_preflight
#' @param ... Ignored.
#' @export
print.mechanismPreflight <- function(x, ...) {
    cat("<mechanism preflight>\n")
    cat("  tier:    ", x$tier, "\n", sep = "")
    m <- x$metrics
    cat(sprintf("  inputs:  %d pathways, %s up genes, median %s members, %s cell line(s)\n",
                m$n_pathways,
                if (is.na(m$n_up)) "?" else as.character(m$n_up),
                if (is.finite(m$median_members)) as.character(as.integer(m$median_members)) else "?",
                if (is.na(m$n_cell_lines)) "?" else as.character(m$n_cell_lines)))
    if (length(x$reasons)) {
        cat("  why:     ", paste(x$reasons, collapse = "\n           "), "\n", sep = "")
    } else {
        cat("  why:      nothing blocked attribution\n")
    }
    invisible(x)
}
