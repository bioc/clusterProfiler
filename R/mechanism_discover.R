## Discovery gate: from one differential expression profile, shortlist the transcription factors
## that could have driven it.
##
## Why this layer. The verification layer (`pathway_evidence()`) answers "given a candidate, is the
## evidence enough"; by design it cannot answer "who", because identifying the driver means testing
## all 265 candidates on every dataset, and at 1-2 seconds per GPSAdb call that is 30-50 minutes per
## profile. So there has to be a cheap first pass: shortlist candidates, then verify each one.
##
## The method is not ours, it is `ulm` (decoupler): take one TF's target weight vector (target = +/-1,
## non-target = 0), correlate it with the whole log2FC profile by Pearson, then turn the sample size
## into a t-value. We tried 12 combinations (ensemble, cross-network, viper, ranking by absolute
## value) and none beat it.
##
## Full-universe benchmark and frozen inputs: dev_docs/mechanism-benchmark-4215.R.
## A missing top-20 rank is a ranking miss, not evidence of network absence.
##
## So nothing here "improves" it; it only does the three things `ulm` does not: **drop self-loops,
## give a direction, disclose how much is covered**.

#' Rank candidate transcription factors for a differential expression profile
#'
#' This is the cheap first stage. It ranks transcription factors by correlating
#' each one's signed target weights against the whole differential expression
#' vector, then leaves the decision to [pathway_evidence()], which checks a
#' named candidate against measured perturbation data.
#'
#' The score is the t-value of that correlation. **A positive score means the
#' factor looks active in this experiment, a negative one means repressed.**
#' The ranking is not a verdict. In a frozen exploratory corpus of 17 profiles,
#' all 17 expected drivers were present and scorable with CollecTRI. Using the
#' default signed ranking, expected-label hit counts were 6/17 at rank one,
#' 10/17 in the top five and 13/17 in the top ten. The corpus includes one
#' knockout and repeated studies/drivers; these are not independent validation
#' estimates or proof of unique causal identification. A rank outside the first
#' twenty must not be counted as network absence. Treat the output as a shortlist.
#'
#' In the truncated top-twenty candidate pool, score and measured target count
#' were positively associated, including among non-reference candidates
#' (Spearman correlation about 0.24). This selected-pool observation does not
#' establish a general size bias or show whether annotation depth, correlated
#' targets or biology explains it. `n_targets` and `share` disclose coverage.
#' No target-count weight is applied: exploratory adjustments have not been
#' independently validated.
#'
#' @param signed Named numeric vector of log2 fold changes, names are gene
#'   symbols. The whole vector is used, not just the significant genes: the
#'   correlation is what separates a factor's targets from everything else.
#' @param network Data frame of transcription factor targets with columns
#'   `source` (factor), `target` (gene symbol) and optionally `weight`
#'   (positive for activation targets, negative for repression targets; 1 when
#'   absent). CollecTRI and DoRothEA can be exported in this shape. Identical
#'   repeated edges are deduplicated; conflicting weights for the same
#'   source-target pair are rejected rather than resolved by row order.
#' @param top Integer. How many candidates to return.
#' @param min_targets Integer. Factors with fewer measured targets than this are
#'   skipped: a correlation over three genes is noise.
#' @param rank_by Character. `"score"` (default) ranks by the signed score, so
#'   factors that look active come first; that is the right default for
#'   treatment and ligand experiments, where the driver is switched on.
#'   `"abs"` ranks by magnitude regardless of sign, which is what a knockdown
#'   or knockout profile needs: its driver scores negative and would otherwise
#'   sink to the bottom. On our benchmark 16 of 17 profiles are activation
#'   types, so `"score"` wins there; that is a property of the corpus, not of
#'   the method.
#' @param keep_self_targets Logical. When `FALSE` (the default) a factor that
#'   lists itself as a target loses that edge. Self-targets inflate the score
#'   for any knockdown or knockout, because the factor's own gene is the
#'   strongest signal in the profile; in CollecTRI 62% of the factors that
#'   reach evaluation carry such an edge, and removing them costs about five
#'   percentage points of ranking performance.
#'
#' @return A data frame with one row per candidate: `rank`, `source`, `score`
#'   (the t-value), `r` (the correlation), `direction` (`"activates"` or
#'   `"represses"`), `n_targets`, and `share` (the fraction of all measured
#'   genes that are targets of this factor; the denominator is `length(signed)`,
#'   not the number of genes in the network).
#'
#' @seealso [pathway_evidence()] to test the shortlist against perturbation
#'   data, and [mechanism_preflight()] to decide whether the profile can carry
#'   any claim at all.
#'
#' @examples
#' set.seed(2)
#' genes <- paste0("g", 1:50)
#' signed <- stats::setNames(c(rep(1.5, 5), stats::rnorm(45, 0, 0.2)), genes)
#' network <- data.frame(source = "TF1", target = genes[1:5], weight = 1)
#' mechanism_discover(signed, network)
#'
#' @export
mechanism_discover <- function(signed, network, top = 10L, min_targets = 5L,
                               rank_by = c("score", "abs"),
                               keep_self_targets = FALSE) {
    rank_by <- match.arg(rank_by)
    if (!is.numeric(signed) || is.null(names(signed)))
        stop("`signed` must be a named numeric vector (gene symbol to log2 fold change).",
             call. = FALSE)
    ## Duplicated gene names misalign `w[match()]` (which fills only the first position) with `x`
    ## (which has a value at every duplicate position), so the correlation score is silently wrong.
    ## The package's `finish()` already deduplicates by largest |lfc|; do the same here.
    if (anyDuplicated(names(signed))) {
        n_dup <- sum(duplicated(names(signed)))
        ord <- order(-abs(signed))
        signed <- signed[ord]
        signed <- signed[!duplicated(names(signed))]
        message(sprintf("Dropped %d duplicated gene name(s), keeping the largest |log2FC|.", n_dup))
    }
    net <- as.data.frame(network, stringsAsFactors = FALSE)
    if (!all(c("source", "target") %in% names(net)))
        stop("`network` must have `source` and `target` columns.", call. = FALSE)
    if (!"weight" %in% names(net)) net$weight <- 1
    net$source <- as.character(net$source)
    net$target <- as.character(net$target)
    net$weight <- as.numeric(as.character(net$weight))
    net <- net[!is.na(net$source) & !is.na(net$target) & is.finite(net$weight), ,
               drop = FALSE]

    ## Self-loops: the TF listed among its own targets. Left in, the TF's own log2FC in a
    ## knockdown/knockout experiment would push its score up directly (62% of TFs in CollecTRI
    ## carry such an edge).
    if (!keep_self_targets) {
        self <- net$source == net$target
        if (any(self)) {
            message(sprintf("Dropped %d self-target edges (%d factor(s)). See `keep_self_targets`.",
                            sum(self), length(unique(net$source[self]))))
        }
        net <- net[!self, , drop = FALSE]
    }
    if (!nrow(net))
        stop("No usable edges left in `network`.", call. = FALSE)

    ## Handling NA: `signed` may hold genes whose value could not be computed. Keeping them turns
    ## the whole correlation into NA, and "this gene has no measured fold change" is not "it did not
    ## change", so drop them and say so.
    if (any(!is.finite(signed))) {
        n_na <- sum(!is.finite(signed))
        message(sprintf("Dropped %d gene(s) with non-finite log2 fold change.", n_na))
        signed <- signed[is.finite(signed)]
    }
    if (!length(signed)) stop("`signed` holds no finite values.", call. = FALSE)

    ## Resolve exact duplicates only; conflicting evidence must not depend on row order.
    net <- net[!duplicated(net[, c("source", "target", "weight")]), , drop = FALSE]
    if (anyDuplicated(net[, c("source", "target")]))
        stop("Conflicting weights for the same source-target pair in `network`; resolve them before scoring.",
             call. = FALSE)
    net <- net[net$target %in% names(signed), , drop = FALSE]
    if (!nrow(net))
        stop("None of the network targets are measured in `signed`.", call. = FALSE)

    ## Every TF is correlated over **the same set of genes**, so df stays constant.
    ##
    ## The gene universe is **the whole measured profile**, not just the network's targets. Two
    ## reasons:
    ##   1. Statistically the question is "are these target genes special relative to the whole
    ##      transcriptome"; targets alone leave no contrast (in a single-TF network w is all 1, the
    ##      variance is 0, and the function refuses outright);
    ##   2. A network is always incomplete, and genes no TF lists are still part of the background.
    ## With broad network coverage (CollecTRI, about 4600 targets) this matches what decoupler
    ## computes on the adj feature.
    genes <- names(signed)
    x <- as.numeric(signed[genes])
    xc <- x - mean(x)
    sx <- sqrt(sum(xc^2))
    if (!is.finite(sx) || sx == 0)
        stop("The measured part of `signed` has no variance.", call. = FALSE)

    split_idx <- split(seq_len(nrow(net)), net$source)
    rows <- lapply(names(split_idx), function(tf) {
        idx <- split_idx[[tf]]
        w <- numeric(length(genes))
        w[match(net$target[idx], genes)] <- net$weight[idx]
        n <- sum(w != 0)
        if (n < min_targets) return(NULL)
        wc <- w - mean(w)
        sw <- sqrt(sum(wc^2))
        if (!is.finite(sw) || sw == 0) return(NULL)
        r <- sum(wc * xc) / (sw * sx)
        r <- max(min(r, 1 - 1e-12), -1 + 1e-12)
        df <- length(genes) - 2L
        data.frame(source = tf, r = r,
                   score = r * sqrt(df / ((1 - r) * (1 + r))),
                   direction = if (r > 0) "activates" else "represses",
                   n_targets = n, share = n / length(genes),
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    if (is.null(out) || !nrow(out))
        stop("No factor had enough measured targets. Lower `min_targets` or widen the network.",
             call. = FALSE)
    ## The number of factors actually scored, not the number with edges in the network. The latter
    ## also counts those skipped by `min_targets`, so "tested 1176" looks like all were scored.
    n_scored <- sum(!vapply(rows, is.null, logical(1)))
    out <- out[order(if (rank_by == "abs") -abs(out$score) else -out$score), ,
               drop = FALSE]
    out$rank <- seq_len(nrow(out))
    out <- utils::head(out, max(1L, as.integer(top)[1]))
    out <- out[, c("rank", "source", "score", "r", "direction", "n_targets",
                   "share"), drop = FALSE]
    rownames(out) <- NULL
    ## Attributes are set after subsetting: `[.data.frame` drops them, so setting them earlier
    ## would be wasted.
    attr(out, "n_genes") <- length(genes)
    attr(out, "n_factors_tested") <- n_scored
    attr(out, "n_factors_in_network") <- length(split_idx)
    class(out) <- c("mechanism_discovery", "data.frame")
    out
}

#' @export
print.mechanism_discovery <- function(x, n = 10L, ...) {
    cat("<mechanism discovery>\n")
    cat("  genes used:  ", attr(x, "n_genes"), "\n", sep = "")
    cat("  factors tested:", attr(x, "n_factors_tested"),
        " (of ", attr(x, "n_factors_in_network"), " in network)\n", sep = "")
    show <- utils::head(as.data.frame(x), n)
    print(show[, c("rank", "source", "score", "direction", "n_targets")],
          row.names = FALSE)
    cat("\nThis is a shortlist, not a verdict: verify the candidates with pathway_evidence().\n")
    invisible(x)
}
