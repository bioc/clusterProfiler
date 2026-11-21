## Pathway-level mechanistic evidence.
##
## Parallel to the TF mode of `evidence()`, but it asks a different question:
##   not "which TF regulates this gene set", but "can this pathway's movement be explained by a
##   measured perturbation".
##
## Designed on `dev_docs/pathway-mechanism-design.md`. Two core constraints:
##
## 1. **Candidate drivers come from the caller** (`drivers`). The tool does not pretend to discover
##    a driver among thousands of datasets -- phase 0 already measured that, and a "who looks most
##    like it" ranking gives no usable answer (biotrigger buried the answer at rank 18/200).
##    What happens here is **a test**: given a candidate, answer "would knocking it down move this
##    pathway".
## 2. **Every conclusion must be refutable**, so the output carries a `statement` line. Something
##    that cannot be falsified is not a mechanistic claim.

## Pathway members and direction.
##
## Direction priority: an explicit `pathway_direction` > inference from the members' log2FC.
## With neither, direction is NA -- **no guessing**. A pathway without direction cannot be tested and
## can only be given level C.
.pathway_members <- function(result) {
    if (is.null(result$geneID)) return(character(0))
    tok <- trimws(unlist(strsplit(as.character(result$geneID), "/")))
    ## Empty tokens and NA are not members. Unfiltered:
    ##   * `"A/B/C/D/E/"` counts 6 members, so `n_members` and the number in the statement are both
    ##     one too high (measured: a trailing slash and a space gave 6 instead of 5);
    ##   * `geneID = NA` yields a member called "NA", so "no members" is reported as 1.
    ## Both move the boundary of the `n_members < 5` `too-small` criterion.
    unique(tok[!is.na(tok) & nzchar(tok)])
}

## Where direction comes from decides whether `support` can be used as evidence.
##   "given"   -- the result states it itself (external information, support is independent)
##   "nes"     -- the sign of GSEA's NES (external information, support is independent)
##   "members" -- inferred from the member majority (**so support is vacuously true**; see
##               .pathway_direction_source)
.pathway_direction_source <- function(result, signed) {
    if (!is.null(result$pathway_direction) && !is.na(result$pathway_direction)) return("given")
    if (!is.null(result$NES) && !is.na(result$NES)) return("nes")
    if (!is.null(signed) && length(signed)) return("members")
    "none"
}

.pathway_direction <- function(result, members, signed) {
    if (!is.null(result$pathway_direction) && !is.na(result$pathway_direction)) {
        return(as.character(result$pathway_direction))
    }
    if (!is.null(result$NES) && !is.na(result$NES)) {
        return(if (result$NES > 0) "up" else "down")
    }
    ## An ORA result carries no direction of its own -- it says the pathway is
    ## enriched, not which way it moved. When the caller supplies signed fold
    ## changes the members answer it: whichever way the majority went is the way
    ## the pathway went. Returning NULL here instead made every row untestable,
    ## which is what the first run of this function did.
    if (is.null(signed) || !length(signed)) return(NULL)
    hit <- members[members %in% names(signed)]
    if (!length(hit)) return(NULL)
    moves <- signed[hit]
    n_up <- sum(moves > 0)
    n_down <- sum(moves < 0)
    if (n_up > n_down) "up" else if (n_down > n_up) "down" else NULL
}

## Member direction agreement: the share of members that move the pathway's way and clear the
## threshold.
##
## It exists for the same reason as `support` in the old spec -- a high score driven by a few genes
## cannot count as support for the whole pathway. When only 5 of a 33-member pathway move the same
## way, that is not "the pathway is moving", it is "a few genes are moving".
.pathway_support <- function(members, signed, direction) {
    if (is.null(signed) || !length(signed) || is.null(direction)) return(NA_real_)
    hit <- members[members %in% names(signed)]
    if (!length(hit)) return(NA_real_)
    moves <- signed[hit]
    want <- if (direction == "up") moves > 0 else moves < 0
    mean(want)
}

## Resolve a candidate driver gene into GPSAdb datasets.
##
## One gene may have several datasets. All of them are fetched and tested one by one, and the
## strongest is used -- but which one was used is recorded, because "the same gene behaves
## differently in different cell lines" really happens (in phase 0, VHL did not match expectation in
## HK-2).
.pathway_datasets <- function(gene, limit = 5L) {
    .pathway_dataset_lookup(gene, limit)$data
}

## Resolve one candidate and say which of three things happened: data was found,
## the lookup succeeded but returned nothing usable, or the lookup itself failed.
##
## The third case used to be indistinguishable from the second. `pathway_evidence()`
## then printed "No usable perturbation data found in GPSAdb" for a network outage
## or a machine without the suggested `gpsadb` package -- a false statement about
## the database, and the same confusion the inventory lookup was fixed for.
.pathway_dataset_lookup <- function(gene, limit = 5L) {
    ## Normalise through the alias table the package already ships, then query.
    ##
    ## Without normalisation the lookup **silently comes back empty**: in GPSAdb `GR` has 0
    ## knockdown datasets and `NR3C1` has 1; a caller writing `drivers = "GR"` is being natural and
    ## gets "this driver has no data at all", which looks exactly like "the data really does not
    ## support it". Measured on airway this is precisely what happens: all 9 pathways C, with no
    ## sign that aliasing was the problem.
    ## The alias table is **one-way** (`NR3C1 -> GR`), while caller and database may each sit at
    ## either end: in GPSAdb `NR3C1` has data and `GR` does not, and users most naturally write
    ## `drivers = "GR"`. So try both directions and use whichever hits first.
    aliases <- tryCatch(mechanism_regulator_aliases(), error = function(e) character(0))
    back <- if (length(aliases)) names(aliases)[aliases == gene] else character(0)
    fwd <- tryCatch(mechanism_alias_regulator(gene), error = function(e) gene)
    tries <- unique(c(gene, fwd, back))
    tries <- tries[!is.na(tries) & nzchar(tries)]

    errors <- character(0)
    for (g in tries) {
        got <- .pathway_lookup_one(g, limit)
        errors <- c(errors, got$errors)
        if (!is.null(got$data)) {
            return(list(data = got$data, status = "available",
                        message = paste(unique(errors), collapse = "; ")))
        }
    }
    list(data = NULL, status = if (length(errors)) "failed" else "empty",
         message = paste(unique(errors), collapse = "; "))
}

## One spelling of the candidate. "knockdown" first because it is the common
## spelling and the narrow query is fast and stable; the other removal spellings
## are only asked when it yields nothing, so a factor recorded solely as CRISPRi
## or knockout (STAT1) is not reported as having no data. An unfiltered query is
## avoided: it is slow and intermittently returns nothing for the largest factors.
.pathway_lookup_one <- function(gene, limit) {
    methods <- c("knockdown", "CRISPRi", "knockout", "siRNA", "shRNA")
    errors <- character(0)
    for (m in methods) {
        got <- tryCatch(
            list(d = as.data.frame(gpsadb::gpsadb_search_datasets(
                perturbed_gene = gene, method = m, limit = 50)[[1]]), error = NULL),
            error = function(e) list(d = NULL, error = conditionMessage(e)))
        if (!is.null(got$error)) {
            errors <- c(errors, sprintf("%s: %s", m, got$error))
            next
        }
        d <- .pathway_usable_datasets(got$d)
        if (is.null(d)) next
        ## Keep one dataset per cell line first (the one with more samples), so the few that are
        ## taken come from different lines -- the reproduction call counts cell lines, not datasets.
        if ("n_samples" %in% names(d)) d <- d[order(-d$n_samples), , drop = FALSE]
        d <- d[!duplicated(d$cell_line), , drop = FALSE]
        return(list(data = utils::head(d, limit), errors = errors))
    }
    list(data = NULL, errors = errors)
}

## Post-query normalisation, factored out so it can be tested offline.
##
## The method filter used to be `method = "knockdown"`. GPSAdb spells removal
## four ways (siRNA, shRNA, CRISPRi, knockout), so that whitelist reported STAT1
## -- which has one CRISPRi and one knockout dataset and no literal "knockdown"
## -- as having no data at all, the same false negative the inventory lookup was
## fixed for. Overexpression is excluded instead.
.pathway_usable_datasets <- function(d) {
    if (is.null(d) || !nrow(d)) return(NULL)
    if (!"cell_line" %in% names(d) || !"dataset_id" %in% names(d)) return(NULL)
    if ("method" %in% names(d)) {
        d <- d[is.na(d$method) | !grepl("overexpress", d$method, ignore.case = TRUE), ,
               drop = FALSE]
    }
    if ("in_reference" %in% names(d)) {
        d <- d[!is.na(d$in_reference) & d$in_reference, , drop = FALSE]
    }
    if (!nrow(d)) return(NULL)
    d$cell_line <- trimws(as.character(d$cell_line))
    d$cell_line[!is.na(d$cell_line) & !nzchar(d$cell_line)] <- NA_character_
    identities <- split(d$cell_line, as.character(d$dataset_id))
    conflicting <- names(identities)[vapply(identities, function(z)
        length(unique(z)) > 1L, TRUE)]
    if (length(conflicting)) {
        warning(sprintf("Ignoring GPSAdb dataset(s) with conflicting cell-line identities: %s.",
                        paste(conflicting, collapse = ", ")), call. = FALSE)
        d <- d[!d$dataset_id %in% conflicting, , drop = FALSE]
    }
    if (!nrow(d)) return(NULL)
    d <- d[!duplicated(d$dataset_id), , drop = FALSE]
    if (!nrow(d)) NULL else d
}

## Pick up the reproduction row on the single-cell-line branch. **Aliases must be resolved on both
## sides**: GPSAdb returns `NR3C1`, while the caller's `mechanism_reproduction()` table may have been
## built with the alias `GR`. Unresolved, the match silently returns 0 rows and the wording quietly
## falls back to the default -- GR happens to be the driver of the two level-B cases in P4, and that
## is where this hole was found.
.pe_reproduction_row <- function(gene, pathway, reproduction) {
    if (is.null(reproduction) || !nrow(reproduction)) return(NULL)
    if (!all(c("driver", "pathway") %in% names(reproduction)))
        stop("`reproduction` must have `driver` and `pathway` columns.", call. = FALSE)
    aliases <- .mechanism_driver_names(gene)
    reproduction[reproduction$driver %in% aliases &
                 as.character(reproduction$pathway) == pathway, , drop = FALSE]
}

## Wording for the single-cell-line branch. Factored out so it **can be tested on its own**: the
## three branches (reproduced / not reproduced / no second dataset) must stay distinguishable, and
## inline inside `pathway_evidence()` they could only be checked by running real data.
##
## Grades come from `mechanism_reproduction()` (spec Section 4.3):
##   reproduces      -> say how many independent datasets reproduced it
##   consistent      -> keep the original wording
##   not-reproduced  -> state plainly that it is unverified
##   not-computable  -> same reading as not-reproduced; must not count as a pass by default
.pe_single_line_statement <- function(gene, n_members, dataset, padj, row) {
    base <- sprintf(
        "Knocking down candidate %s moves these %d members opposite the observed pathway direction in %s (padj %.2g), with one supporting cell line.",
        gene, n_members, dataset, padj)
    grade <- if (is.null(row) || !nrow(row)) NA_character_ else as.character(row$grade[1])
    if (!is.na(grade) && grade == "reproduces") {
        agree <- as.integer(round(row$concordance[1] * row$k[1]))
        sprintf("%s It also reproduces across %d independent datasets (%d/%d agree).",
                base, as.integer(row$n_biosamples[1]), agree, as.integer(row$k[1]))
    } else if (!is.na(grade) && grade == "not-reproduced") {
        agree <- as.integer(round(row$concordance[1] * row$k[1]))
        sprintf("%s The available perturbation data do not reproduce that shift (%d/%d agree), so treat this as unverified.",
                base, agree, as.integer(row$k[1]))
    } else if (!is.na(grade) && grade == "not-computable") {
        sprintf("%s No independent reproduction is available; this is consistent perturbation evidence, but under-evidenced, not identification of a driver.",
                base)
    } else {
        sprintf("%s This is a limit of the available data: consistent perturbation evidence, but under-evidenced, not identification of a driver.",
                base)
    }
}

#' Test enriched pathways against measured perturbation data
#'
#' @param x An `enrichResult` or `gseaResult` from [enrichGO()], [enrichKEGG()]
#'   and friends, or a `$pathways` table from an [evidence()] bundle.
#' @param drivers Character vector of candidate driver genes. Each is resolved to
#'   its perturbation datasets and tested against every pathway.
#' @param signed Optional named numeric of log2 fold changes, used to work out
#'   which way a pathway's members moved when the result does not say.
#' @param alpha Significance threshold for a perturbation to count.
#' @param min_support Minimum fraction of members moving the same way before a
#'   pathway counts as moved at all.
#' @param max_pathways How many pathways to report, taken in the result's own order.
#'   This also caps the run: every pathway is tested against every driver dataset
#'   over the network, so the wall time grows with `max_pathways * drivers *
#'   datasets`. Measured on the hypoxia pilot, 25 pathways took 96 calls and about
#'   a minute, 60 took 210 calls and about three. Results are memoised by
#'   (dataset, member set), which helps when pathways share members and not
#'   otherwise, so pick this deliberately rather than leaving it high.
#' @param reproduction Optional data frame from [mechanism_reproduction()]. When
#'   given, a level-B row whose only support is a single cell line is worded
#'   according to how many independent perturbation datasets agree: a pair that
#'   reproduces says so, a pair that does not reproduce is labelled unverified,
#'   and a pair with no second dataset remains consistent but under-evidenced.
#' @param max_datasets_per_driver Cap on how many datasets are tested per driver.
#'   Testing all of them is `n_drivers * n_datasets` calls and grows without
#'   bound; six is enough to see two cell lines agree, which is all level A
#'   needs. Datasets from distinct cell lines are preferred over extra replicates
#'   of one line, because replication is what the grade turns on. The cap applies
#'   to datasets tested, not to the lookups that resolve them: a candidate is
#'   searched under one or more removal spellings, so an unknown symbol can cost
#'   a handful of search calls before it is reported as empty.
#' @param n For `print()`, how many rows to show.
#' @param ... Ignored, for method consistency.
#' @return A `data.frame`, one row per pathway, classed `pathway_evidence` so it
#'   prints as a short summary rather than twenty columns. Two columns carry the
#'   interpretation: `statement` holds a claim a perturbation experiment could
#'   refute, and `reason` says why a row stopped short of level A
#'   (`single-cell-line`, `direction-conflict`, `no-dataset`, `too-small`,
#'   `not-significant`, `member-derived` and the input-side cases). Level A leaves
#'   `reason` as `NA`. Level A is weaker than its name suggests: a perturbation that
#'   moves the pathway also moves it if the perturbation is unrelated but perturbs
#'   transcription broadly, and four of thirty background drivers reached A on the
#'   hypoxia pilot. Read A as consistency worth following up, not as identification
#'   of the driver; the statement says as much. A level-B row distinguishes `coherent-untested` (no test
#'   could run) from the legacy key `coherent-refuted` (tests ran but none reached
#'   significance; this is not evidence of falsity). When the answer cannot reach A
#'   for a structural reason the call also warns.
#'   The `perturbation_tests` attribute is a data frame with character columns
#'   `pathway`, `driver`, `dataset`, `cell_line`, `status`, `direction` and double
#'   columns `es`, `padj`. It retains each candidate/dataset outcome, including
#'   `supports` (significant, opposite the pathway direction), `conflicts`
#'   (significant, same direction), `not-significant`, `failed`, and
#'   `not-computable` (missing direction, missing/nonfinite statistics or zero ES).
#'   These statuses describe consistency, not unique causal identification.
#'   A known member support at or above `min_support` is required for promotion
#'   based on significant perturbation evidence. Replication is counted within
#'   the selected candidate only; contradictory datasets remain in the attribute.
#'   The `dataset_lookup` attribute records, per requested driver, whether the
#'   lookup returned data (`available`), succeeded but found nothing usable
#'   (`empty`), or could not be performed (`failed`, with the error text). A
#'   failed lookup yields reason `lookup-failed` and is never reported as the
#'   database holding no data; `not-run` marks an empty input that never reached
#'   dataset resolution. With several drivers, a partial failure is disclosed in
#'   that attribute and by a warning rather than in the row `reason`, which keeps
#'   describing the candidate that was actually tested.
#' @export
pathway_evidence <- function(x, drivers = NULL, signed = NULL, alpha = 0.05,
                             min_support = 0.5, max_pathways = 20L,
                             max_datasets_per_driver = 6L,
                             reproduction = NULL) {
    d <- if (is.data.frame(x)) x else as.data.frame(x)
    if (!nrow(d)) {
        ## An empty input must still speak up. Without this line the user gets an empty data.frame
        ## and cannot tell "the enrichment genuinely had no result" from "the function broke" --
        ## measured, two datasets (estradiol, and other condition experiments that enrich nothing)
        ## stop right here.
        warning(paste("The enrichment result holds no pathways, so there is nothing to score.",
                      "That is usually the condition experiment failing to concentrate its response,",
                      "not a fault in pathway_evidence."), call. = FALSE)
        return(.pathway_evidence_empty(drivers))
    }
    ## The guard sits **after** `head()`. Placed before it, `max_pathways = 0L` left d with 0 rows,
    ## the later `do.call(rbind, list())` returned NULL, and `attr(out, ...)` failed with
    ## "attempt to set an attribute on NULL" -- an error that names no cause.
    max_pathways <- max(1L, as.integer(max_pathways)[1])
    d <- utils::head(d, max_pathways)
    if (!nrow(d)) return(.pathway_evidence_empty(drivers))

    ## Identifier space: **detect and translate automatically**, instead of making the caller bitr
    ## every time.
    ##
    ## ReactomePA returns Entrez, GPSAdb only understands symbols, and the combination of those two
    ## sources is the main path. Earlier the caller had to convert, at the cost of anyone who forgot
    ## getting a "ran fine, nothing happened" result -- 0.8 seconds, zero network calls, everything
    ## untestable. Detection is cheap (do the members look purely numeric), so doing it
    ## automatically is more reliable than asking people to remember.
    d <- .pathway_ids_to_symbols(d, signed)

    ## Identifier-space guard.
    ##
    ## ReactomePA's `geneID` is **Entrez**, while GPSAdb and most differential tables use **Symbols**.
    ## When the two do not line up, no member of any pathway can match `signed`, direction cannot be
    ## inferred, support is NA, and the result is **everything untestable while looking like a normal
    ## null conclusion** -- measured, it finished in 0.8 seconds without sending a single network
    ## call. That "ran fine, nothing happened" is far more dangerous than an error.
    ##
    ## The criterion is "not one pathway matched": as soon as one matches, the spaces are connected,
    ## and a few unmatched pathways are normal (their members are not in the differential table).
    if (!is.null(signed) && length(signed) && nrow(d)) {
        any_hit <- any(vapply(seq_len(nrow(d)), function(i) {
            m <- .pathway_members(d[i, ])
            length(m) && any(m %in% names(signed))
        }, TRUE))
        if (!any_hit) {
            rlang::abort(c(
                "No pathway member matched the `signed` gene list.",
                x = "Every pathway came back untestable, which usually means the two are in different identifier spaces.",
                i = "ReactomePA returns Entrez IDs; pass symbols in both, or translate one side first."
            ))
        }
    }

    ## Candidate datasets are resolved once, not re-queried per pathway -- one driver has the same
    ## datasets across dozens of pathways.
    ## `lookup` records each candidate's resolution: data found, succeeded but empty, or the query
    ## itself failed. Conflating the three would print "GPSAdb has no data for this factor" after a
    ## network outage.
    sets <- list()
    lookup <- list()
    if (length(drivers)) {
        for (g in drivers) {
            lk <- .pathway_dataset_lookup(g, limit = max_datasets_per_driver)
            lookup[[g]] <- data.frame(driver = g, status = lk$status,
                message = substr(lk$message, 1, 300), stringsAsFactors = FALSE)
            if (!is.null(lk$data))
                sets[[g]] <- lk$data[, c("dataset_id", "cell_line"), drop = FALSE]
        }
    }
    lookup <- if (length(lookup)) do.call(rbind, lookup) else
        data.frame(driver = character(), status = character(), message = character(),
                   stringsAsFactors = FALSE)
    failed <- lookup$driver[lookup$status %in% "failed"]

    ## That level A is unreachable can be settled **at the driver level**, without back-inferring it
    ## from each pathway's `reason`.
    ##
    ## The earlier version waited for the results and then looked at the distribution of reasons,
    ## and it under-reported: half the airway pathways were direction conflicts and half were split
    ## members, so there was no "single dominant reason" and GR's single cell line went unreported.
    ## The criterion belongs here, where it is both more direct and not prone to false alarms.
    if (length(sets)) {
        cells_each <- vapply(sets, function(z) {
            cells <- trimws(as.character(z$cell_line))
            length(unique(cells[!is.na(cells) & nzchar(cells)]))
        }, 0L)
        if (all(cells_each < 2L)) {
            warning(sprintf(
                "Level A is unreachable: %s has perturbation data in a single cell line, and level A requires two. Results here cannot exceed B.",
                paste(sprintf("`%s` (%d)", names(cells_each), cells_each), collapse = ", ")),
                call. = FALSE)
        }
    } else if (length(drivers)) {
        if (length(failed)) {
            warning(sprintf(paste(
                "GPSAdb could not be queried for %s; no dataset was resolved.",
                "This is a lookup failure (network or missing gpsadb), not evidence that these factors have no data.",
                "Every pathway will be described without perturbation evidence."),
                paste0("`", paste(failed, collapse = "`, `"), "`")), call. = FALSE)
        } else {
            warning(sprintf(
                "No usable perturbation data found in GPSAdb for %s; every pathway will be described without perturbation evidence.",
                paste0("`", paste(drivers, collapse = "`, `"), "`")), call. = FALSE)
        }
    }
    if (length(sets) && length(failed)) {
        warning(sprintf(
            "GPSAdb could not be queried for %s; those candidates are not tested. This is a lookup failure, not evidence that they have no data.",
            paste0("`", paste(failed, collapse = "`, `"), "`")), call. = FALSE)
    }

    ## Memoise by (dataset, member signature).
    ##
    ## Caching per pathway is wrong -- the result depends on the members, and a different pathway is
    ## a different test. But **different pathways often share exactly the same member set**: in the
    ## hypoxia result `response to hypoxia`, `response to decreased oxygen levels` and `response to
    ## oxygen levels` are the same 33 genes and used to be tested three times over. Reusing on an
    ## identical signature is safe because the request is identical.
    memo <- new.env(parent = emptyenv())
    calls <- 0L
    fails <- 0L
    perturbation_tests <- list()
    record_test <- function(pathway, gene, id, cell, direction, got = NULL) {
        es <- if (is.null(got$es)) NA_real_ else got$es
        padj <- if (is.null(got$padj)) NA_real_ else got$padj
        status <- if (isTRUE(got$.failed)) "failed"
            else if (is.null(direction) || !is.finite(es) || !is.finite(padj) || es == 0) "not-computable"
            else if (padj >= alpha) "not-significant"
            else if ((direction == "up" && es < 0) ||
                     (direction == "down" && es > 0)) "supports" else "conflicts"
        perturbation_tests[[length(perturbation_tests) + 1L]] <<- data.frame(
            pathway = pathway, driver = gene, dataset = id, cell_line = cell,
            es = as.double(es), padj = as.double(padj), status = status,
            direction = if (is.null(direction)) NA_character_ else direction,
            stringsAsFactors = FALSE)
    }
    ask <- function(id, members) {
        key <- paste0(id, "|", paste(sort(members), collapse = ","))
        if (!is.null(memo[[key]])) return(memo[[key]])
        calls <<- calls + 1L
        val <- tryCatch({
            r <- gpsadb::gpsadb_test_gene_set(id, members)
            t <- if (length(r)) as.data.frame(r[[1]]) else data.frame()
            if (!nrow(t)) list(padj = NA_real_, es = NA_real_) else
                list(padj = if (length(t$padj)) as.numeric(t$padj[1]) else NA_real_,
                     es = if (length(t$es)) as.numeric(t$es[1]) else NA_real_)
        }, error = function(e) list(.failed = TRUE))
        if (isTRUE(val$.failed)) fails <<- fails + 1L
        memo[[key]] <- val
        val
    }

    rows <- lapply(seq_len(nrow(d)), function(i) {
        members <- .pathway_members(d[i, ])
        dsrc <- .pathway_direction_source(d[i, ], signed)
        direction <- .pathway_direction(d[i, ], members, signed)
        support <- .pathway_support(members, signed, direction)
        desc <- .pathway_label(d[i, ])
        ## Default first: the "no direction" branch below returns early, before any test has run,
        ## yet it still has to fill the `n_tests_run` column.
        n_tested <- 0L

        ## No direction has two causes, and they must be two conclusions, not one:
        ##   * `signed` was not supplied -- we never measured which way the members moved, this is
        ##     **untested**;
        ##   * `signed` was supplied but the members split evenly -- measured, they have no common
        ##     direction, this is **tested and not supported**.
        ## Merging them into one "no direction" calls the latter missing, and keeping it apart from
        ## `absent` / `untested` is the same discipline.
        if (is.null(direction)) {
            for (g in names(sets)) {
                for (j in seq_len(nrow(sets[[g]]))) {
                    record_test(desc, g, as.character(sets[[g]]$dataset_id[j]),
                                as.character(sets[[g]]$cell_line[j]), direction)
                }
            }
            untested <- is.null(signed) || !length(signed)
            ## Here there are really **three** conclusions, not one:
            ##   * no `signed`            -> no-fold-changes (untested)
            ##   * supplied, no member in it -> no-members (these genes are not in the differential
            ##     table)
            ##   * supplied, members present, split evenly -> no-direction (tested, no common
            ##     direction)
            ## The last two used to share "Members split evenly", so a pathway with 0 members (e.g.
            ## `geneID` NA or all unrecognised names) was also called "members split evenly" -- the
            ## same class of error as folding member-derived into low-support in `why_c`.
            hit <- if (untested) character(0) else members[members %in% names(signed)]
            why <- if (untested) "no-fold-changes"
                else if (!length(hit)) "no-members"
                else "no-direction"
            return(data.frame(
                pathway = desc, direction = NA_character_, n_members = length(members),
                support = support, driver = NA_character_, dataset = NA_character_,
                cell_line = NA_character_, es = NA_real_, padj = NA_real_,
                effect = NA_character_, n_tests_run = n_tested,
                n_cell_lines_supported = 0L, level = "C",
                verdict = if (untested) "untestable" else "unsupported",
                reason = why,
                statement = switch(why,
                    "no-fold-changes" = "No signed fold changes were supplied, so the pathway's movement is unknown.",
                    "no-members" = "None of this pathway's members appear in the signed gene list, so its movement is unknown.",
                    "no-direction" = "Members split evenly, so the pathway has no common direction to explain."),
                stringsAsFactors = FALSE))
        }

        ## Test **every** dataset of every candidate driver, not just the strongest one.
        ##
        ## The original kept only the smallest padj, so level A meant "one dataset supports it".
        ## Measured, that meaning is far too weak: three of four drivers stopped holding in another
        ## cell line, and the effect sign flipped (`results/stage3_refutation.tsv`). All datasets
        ## have to be visible for reproduction to be judgeable.
        tested <- list()
        n_attempt <- 0L
        n_fail <- 0L
        for (g in names(sets)) {
            for (j in seq_len(nrow(sets[[g]]))) {
                id <- as.character(sets[[g]]$dataset_id[j])
                got <- ask(id, members)
                record_test(desc, g, id, as.character(sets[[g]]$cell_line[j]), direction, got)
                n_attempt <- n_attempt + 1L
                ## A failed call and a call that succeeded with no result must stay apart.
                ## Both used to be `next`, so one outage emptied `tested`, which was then read as
                ## "the members disagree in direction" -- measured, HIF1A (4 cell lines) came back
                ## as `low-support` under a simulated outage with zero warnings, and the user could
                ## not see it was the network.
                if (isTRUE(got$.failed)) { n_fail <- n_fail + 1L; next }
                if (!is.finite(got$padj) || !is.finite(got$es) || got$es == 0) next
                ## The same dataset may be reached through two drivers aliasing the same gene; record
                ## it once.
                already <- vapply(tested, function(z) identical(z$dataset, id) &&
                                  identical(z$gene, g), TRUE)
                if (any(already)) next
                tested[[length(tested) + 1L]] <- list(
                    gene = g, dataset = id,
                    cell_line = as.character(sets[[g]]$cell_line[j]),
                    es = got$es, padj = got$padj)
            }
        }
        ## The strongest significant result is the one a claim is made on; reproduction requires
        ## that **another cell line** is also significant and in the same direction.
        best <- NULL
        for (cand in tested) {
            if (is.na(cand$padj) || cand$padj >= alpha) next
            if (is.null(best) || cand$padj < best$padj) best <- cand
        }
        n_tested <- length(tested)

        ## A test that could not run has two external causes, and both must **override**
        ## `no-perturbation-test`:
        ##   * every call failed (outage, upstream down) -- a temporary fault, not a property of the
        ##     data;
        ##   * the driver has no dataset at all in GPSAdb -- a data supply problem the user can fix.
        ## Both used to be caught first by the B branch on pathways with external direction and
        ## enough support, and all of them were called "no perturbation test identified a driver" --
        ## measured, that is what an outage printed, with no sign it was the network.
        blocked <- if (n_attempt > 0L && n_fail == n_attempt) "no-test-completed"
            else if (length(drivers) && !length(sets))
                if (length(failed)) "lookup-failed" else "no-dataset"
            else NA_character_

        cells_agree <- character(0)
        if (!is.null(best)) {
            want <- sign(best$es)
            cells_agree <- unique(vapply(tested, function(cand) {
                if (!identical(cand$gene, best$gene)) return(NA_character_)
                if (is.na(cand$padj) || cand$padj >= alpha) return(NA_character_)
                if (is.na(cand$es) || sign(cand$es) != want) return(NA_character_)
                cand$cell_line
            }, character(1)))
            cells_agree <- trimws(cells_agree)
            cells_agree <- unique(cells_agree[!is.na(cells_agree) & nzchar(cells_agree)])
        }
        replicated <- length(cells_agree) >= 2L

        if (!is.null(best) && best$padj < alpha &&
            (!is.na(support) && support >= min_support)) {
            ## Significance is not enough: **the direction must line up**.
            ## If the pathway runs up and its driver is knocked down, the pathway should **fall
            ## back** (es < 0) -- that says the factor switches it on. If knocking it down pushes the
            ## pathway **further up** (es > 0) then the factor was repressing it, not driving it.
            ## The original ignored this and recorded any significant result as support.
            effect <- if (is.na(best$es) || is.na(direction)) NA_character_ else
                if ((direction == "up" && best$es < 0) ||
                    (direction == "down" && best$es > 0)) "activator" else "repressor"
            agrees <- is.na(effect) || effect == "activator"
            ## Level A demands **reproduction**: at least two different cell lines, both significant
            ## and in the same direction.
            ## Single-dataset support drops to B and is marked `single-dataset` -- not because there
            ## is no evidence, but because the evidence has one source, and measured, 3/4 of such
            ## evidence flips when the data changes.
            level <- if (agrees && replicated) "A" else "B"
            verdict <- if (!agrees) "direction-conflict"
                       else if (replicated) "perturbation-supported"
                       else "single-dataset"
            ## A level-A statement must name what it did **not** establish. Measured, on the same
            ## hypoxia pathway, 4 of 30 drivers drawn at random from the inventory (all essential
            ## genes of the splicing / RNA processing / Polycomb kind) also reached A -- knocking down
            ## an essential gene derails the transcriptome and any pathway of a dozen to a few dozen
            ## members moves with it. So "knocking it down moves this" is true, "it is the driver"
            ## is not.
            statement <- if (agrees && replicated) sprintf(
                "If %s drives this programme, knocking it down should move these %d members; %s does (padj %.2g) and %d cell lines agree. Consistent with drive, but not evidence that %s is the driver.",
                best$gene, length(members), best$dataset, best$padj, length(cells_agree), best$gene)
            ## With single-dataset support the reading must be spelled out: **this is not enough
            ## data, not evidence against**. Measured (11 known-truth datasets): 2 of the level-B rows
            ## pointed at the correct driver (GR in airway and GSE286890), and both rows held down by
            ## the "single cell line" demotion were right. The old wording was only a demotion note,
            ## so readers took it as defective.
            else if (agrees) .pe_single_line_statement(
                best$gene, length(members), best$dataset, best$padj,
                .pe_reproduction_row(best$gene, desc, reproduction))
            else sprintf(
                "%s moves these %d members in the same direction the pathway already runs (padj %.2g), consistent with repression rather than drive, not identification of a causal role.",
                best$gene, length(members), best$padj)
            conflicts <- sum(vapply(tested, function(cand) {
                identical(cand$gene, best$gene) && cand$padj < alpha &&
                    sign(cand$es) != sign(best$es)
            }, TRUE))
            if (conflicts > 0L) statement <- paste0(statement, sprintf(
                " However, %d significant dataset(s) for this same candidate have the opposite effect; the evidence is contradictory across datasets.", conflicts))
            return(data.frame(
                pathway = desc, direction = direction, n_members = length(members),
                support = support, driver = best$gene, dataset = best$dataset,
                cell_line = best$cell_line, es = best$es, padj = best$padj,
                effect = effect, n_tests_run = n_tested,
                n_cell_lines_supported = length(cells_agree),
                level = level, verdict = verdict,
                ## The demotion reason must be **written out**. `single-dataset` and
                ## `direction-conflict` both land on B, but one is "not enough evidence" and the
                ## other "the evidence points the opposite way" -- the reader has to be able to
                ## separate them.
                reason = if (level == "A") NA_character_
                    else if (!agrees) "direction-conflict" else "single-cell-line",
                statement = statement,
                stringsAsFactors = FALSE))
        }

        ## Without measured support, agreeing members are only weak evidence (level B) -- and only
        ## **when the direction comes from outside**.
        ##
        ## If direction is inferred from the member majority then `support` is vacuously true (a
        ## majority is necessarily >= 50%), so using it as evidence both fixes the direction and
        ## proves it with the same genes. Found in review: `min_support = 0.5` always held, which left
        ## level B with no information.
        if (is.na(blocked) && !is.na(support) && support >= min_support && dsrc != "members") {
            ## A test that **ran** and found nothing significant, and one that **never ran**, differ.
            ## Both used to share `coherent-untested` / `no-perturbation-test`, so a pathway where 4
            ## datasets were tested and none was significant reported "no perturbation test could
            ## point at a driver" -- the reader could not tell "tested and not supported" from
            ## "never tested", and the word `untested` is simply wrong once testing has happened.
            ran <- length(tested) > 0L
            return(data.frame(
                pathway = desc, direction = direction, n_members = length(members),
                support = support, driver = NA_character_, dataset = NA_character_,
                cell_line = NA_character_, es = NA_real_, padj = NA_real_,
                effect = NA_character_, n_tests_run = n_tested,
                n_cell_lines_supported = 0L, level = "B",
                verdict = if (ran) "coherent-refuted" else "coherent-untested",
                reason = if (ran) "not-significant" else "no-perturbation-test",
                statement = if (ran) sprintf(
                    "%.0f%% of members moved %s, so the pathway moved; the %d perturbation test(s) run against it found no significant movement.",
                    100 * support, direction, n_tested)
                else sprintf(
                    "%.0f%% of members moved %s, so the pathway moved; no perturbation test could be run.",
                    100 * support, direction),
                stringsAsFactors = FALSE))
        }

        ## A C row that gets here has three completely different meanings, and they must be told
        ## apart:
        ##   * members are not in the differential table          -> no-members
        ##   * members do not line up, direction has no external basis -> low-support / member-derived
        ##   * tests ran but none produced a p-value             -> too-small (pathway too small)
        ##   * tests ran and none was significant                -> not-significant
        ## All of them used to be `C unsupported`, and the reader could not tell "not found" from
        ## "could not be looked for".
        ## The order matters: when members number in the single digits, **that is the reason to
        ## report**. A 2-member pathway can only take its direction from those 2 genes, so it
        ## necessarily lands on `low-support` -- but a reader seeing "members disagree" will go and
        ## change the direction strategy when the real problem is too few differential genes.
        ## Measured, MCF-7 (23 up-regulated) and A549 etoposide (16 up-regulated) were both
        ## mislabelled low-support for this reason.
        ## `member-derived` and `low-support` **must stay apart**. They were merged into one
        ## "low-support", so a pathway with support = 100% printed "below the 50% floor" -- I thought
        ## that was fixed last round, but the `else` branch added for it could never be reached,
        ## because why_c never took that value.
        ## This time every value has a definite source and `statement` maps one-to-one through
        ## `switch`, with no catch-all.
        ## `no-test-completed` is ordered before `member-derived`: when every call collapses, "the
        ## test could not run" deserves to be seen more than "direction was inferred from members" --
        ## otherwise users think their input is wrong when it is the network. Measured during an
        ## outage, the per-row reason once showed `member-derived` and only the global warning told
        ## the truth, so the two layers disagreed.
        ## `support` is NA for two different inputs: no `signed` list at all, or
        ## members that are not in it. Reporting both as "members were not found
        ## in the signed gene list" misdiagnosed the first case -- a documented
        ## way to call the function when direction comes from `pathway_direction`
        ## or NES -- and told the user to check identifiers when the missing
        ## input was the fold changes.
        missing_fc <- is.null(signed) || !length(signed)
        why_c <- if (is.na(support)) {
                if (missing_fc) "no-fold-changes" else "no-members"
            }
            else if (!is.na(blocked)) blocked
            else if (length(members) < 5L) "too-small"
            else if (support < min_support) "low-support"
            else if (identical(dsrc, "members")) "member-derived"
            else if (!length(tested)) {
                if (length(sets)) "too-small" else "no-dataset"
            } else "not-significant"

        data.frame(
            pathway = desc, direction = direction, n_members = length(members),
            support = support, driver = NA_character_, dataset = NA_character_,
            cell_line = NA_character_, es = NA_real_, padj = NA_real_,
            effect = NA_character_, n_tests_run = n_tested,
            n_cell_lines_supported = 0L, level = "C", verdict = "unsupported",
            reason = why_c,
            ## The statement must agree with `why_c`. There used to be only "below the support
            ## floor", but support = 100% also comes through here (because the direction was
            ## inferred from the members themselves), so it printed the self-contradiction "100% of
            ## members moved up, below the 50% floor". Measured, that is what airway printed.
            statement = switch(why_c,
                "no-members" = "Members were not found in the signed gene list; the pathway's movement is unknown.",
                "no-fold-changes" = "No signed fold changes were supplied, so member support is unknown; perturbation evidence is retained in the ledger but cannot establish the pathway's direction.",
                "too-small" = sprintf("Only %d members reached this pathway, too few to establish its direction or test a driver; supply more differentially expressed genes.",
                                      length(members)),
                "member-derived" = "The pathway's direction was inferred from these same members, so majority agreement cannot serve as evidence for it.",
                "low-support" = sprintf("%.0f%% of members moved %s, below the support floor of %.0f%%.",
                                        100 * support, direction, 100 * min_support),
                "no-test-completed" = "Every perturbation call for this pathway failed, so the absence of evidence may be an outage rather than a result.",
                "lookup-failed" = "GPSAdb could not be queried for the candidate drivers, so no perturbation test could run; this is a lookup failure, not evidence about the factor.",
                "no-dataset" = "No usable perturbation dataset was available for the candidate drivers, so the pathway was described but not tested.",
                "not-significant" = sprintf("The candidates were tested against these %d members and none reached significance; this does not establish that their effects are absent.",
                                            length(members)),
                "unsupported"),
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, rows)

    ## The cross-source concordance endpoint was removed, for the reason written in the design
    ## document Section 6.4: all three second sources (PPI density, member overlap, the official
    ## GO<->Reactome mapping) failed when measured, and on the main path this column is always
    ## `no-crossref` -- it changes no verdict, and keeping it would make readers think a validation
    ## had been done. **A constant null conclusion misleads more than a missing column.**
    ## Restoring it needs a second source that holds up first, not this column first.

    out <- out[order(match(out$level, c("A", "B", "C")), -out$support), , drop = FALSE]
    attr(out, "drivers") <- drivers
    attr(out, "dataset_lookup") <- lookup
    attr(out, "network_calls") <- calls
    attr(out, "perturbation_tests") <- if (length(perturbation_tests))
        do.call(rbind, perturbation_tests) else .pathway_perturbation_tests_empty()

    ## A structural ceiling must be **said out loud**, not just show up in the grade.
    ##
    ## Without this block, "all B" / "all C" looks exactly like "tested and not supported" -- and
    ## measured, five of nine times the low grade came from externally missing data (GR has only 1
    ## cell line, the pathways are too small, too few differential genes), not from a judgement the
    ## method made. A reader who cannot see the difference reads "cannot be assessed" as "no
    ## evidence".
    ## When every call failed, say so directly without looking at the reason distribution -- that is
    ## a fault at the call level, and no grade should be taken seriously.
    outage <- calls > 0L && fails == calls
    if (outage) {
        warning(sprintf("All %d GPSAdb calls failed. No perturbation evidence was collected; every grade below reflects an outage, not the data. Check the connection and rerun.", calls),
                call. = FALSE)
    } else if (fails > 0L) {
        warning(sprintf("%d of %d GPSAdb calls failed; perturbation evidence is incomplete. Inspect attr(result, 'perturbation_tests') and rerun.", fails, calls),
                call. = FALSE)
    }
    ## No second message during an outage: the `_pathway_warn_capped` line says the same thing.
    ## Two synonymous warnings for one fault only make it look like two problems.
    if (!outage) .pathway_warn_capped(out, names(sets), failed)
    ## The class must be set, otherwise `print.pathway_evidence` is never dispatched.
    ## Measured, this code **never ran once** from the day it was written: a bare data.frame came
    ## back, the user saw the raw 20-column table, and none of the `reason` / `capped by` lines
    ## printed. Another "looks written, never actually ran" example.
    class(out) <- c("pathway_evidence", "data.frame")
    out
}

.pathway_perturbation_tests_empty <- function() {
    data.frame(pathway = character(), driver = character(), dataset = character(),
               cell_line = character(), es = double(), padj = double(),
               status = character(), direction = character(), stringsAsFactors = FALSE)
}

.pathway_evidence_empty <- function(drivers = NULL) {
    out <- data.frame(
        pathway = character(), direction = character(), n_members = integer(),
        support = numeric(), driver = character(), dataset = character(),
        cell_line = character(), es = numeric(), padj = numeric(),
        effect = character(), n_tests_run = integer(),
        n_cell_lines_supported = integer(),
        level = character(), verdict = character(), reason = character(),
        statement = character(), stringsAsFactors = FALSE)
    attr(out, "drivers") <- drivers
    attr(out, "network_calls") <- 0L
    attr(out, "perturbation_tests") <- .pathway_perturbation_tests_empty()
    attr(out, "dataset_lookup") <- data.frame(driver = drivers %||% character(0),
        status = rep("not-run", length(drivers %||% character(0))),
        message = rep("", length(drivers %||% character(0))), stringsAsFactors = FALSE)
    ## Same class as a normal return. Without it, an empty result skips print.pathway_evidence and
    ## the user sees two output shapes from the same function.
    class(out) <- c("pathway_evidence", "data.frame")
    out
}

## A column subset must fall back to a plain data.frame.
##
## `[.data.frame` carries the custom class along, so `out[, c("pathway","level")]` still claims to
## be `pathway_evidence` while `print.pathway_evidence` reads `driver`, `reason` and `support`,
## columns that are no longer there; `is.na(NULL)` is `logical(0)`, and `if` given length 0 fails
## with "argument is of length zero" -- an error that names no cause.
##
## Column subsetting is extremely common (it is everywhere in this project's own scripts), so the
## class's survival is tied to "the columns are all still there" rather than left to print to guess.
#' @export
`[.pathway_evidence` <- function(x, i, j, drop = TRUE) {
    needed <- c("pathway", "direction", "n_members", "support", "driver",
                "cell_line", "es", "padj", "level", "reason", "statement")
    out <- NextMethod()
    if (is.data.frame(out)) {
        attr(out, "drivers") <- attr(x, "drivers", exact = TRUE)
        attr(out, "dataset_lookup") <- attr(x, "dataset_lookup", exact = TRUE)
        attr(out, "network_calls") <- attr(x, "network_calls", exact = TRUE)
        ## A row subset and a column subset lose attributes differently:
        ## `[.data.frame` keeps them for row subsets, but a column subset --
        ## including a pure reorder such as `pe[, rev(names(pe))]` -- drops every
        ## custom attribute while keeping the class. That silently removed the
        ## retained evidence ledger, so `mechanism_report()` on a reordered table
        ## reported `perturbation_tests = NULL` and no candidate evidence.
        ## Whenever the pathway column survives, restore the ledger scoped to the
        ## rows that remain; without it the attribute cannot be scoped and the
        ## class is dropped below.
        ledger <- attr(x, "perturbation_tests", exact = TRUE)
        if (is.data.frame(ledger)) {
            paths <- if ("pathway" %in% names(out)) out$pathway else NULL
            if (is.null(paths)) {
                ## Row-and-column subset that dropped the pathway column: the
                ## kept rows still identify the ledger rows.
                row_subset <- !missing(i) && (nargs() - !missing(drop) >= 3L)
                if (row_subset && !missing(i)) paths <- x$pathway[i]
            }
            if (!is.null(paths))
                attr(out, "perturbation_tests") <-
                    ledger[ledger$pathway %in% paths, , drop = FALSE]
        }
    }
    if (is.data.frame(out) && !all(needed %in% names(out))) {
        class(out) <- "data.frame"
    }
    out
}

## Presentation only: the retained ledger, never A/B/C, determines evidence state.
.pathway_evidence_status <- function(x) {
    tests <- attr(x, "perturbation_tests", exact = TRUE)
    lookup <- attr(x, "dataset_lookup", exact = TRUE)
    valid <- is.data.frame(tests) && all(c("pathway", "status") %in% names(tests))
    empty <- data.frame(evidence_status = character(), availability = character(),
        incomplete = logical(), n_attempted = integer(), n_computable = integer(),
        n_support = integer(), n_conflict = integer(), n_nonsignificant = integer(),
        n_noncomputable = integer(), n_failed = integer(), stringsAsFactors = FALSE)
    if (!nrow(x)) return(empty)
    rows <- lapply(as.character(x$pathway), function(path) {
        s <- if (valid) as.character(tests$status[tests$pathway %in% path]) else character()
        ns <- sum(s %in% "supports"); nc <- sum(s %in% "conflicts")
        nn <- sum(s %in% "not-significant"); nx <- sum(s %in% "not-computable")
        nf <- sum(s %in% "failed")
        ncomp <- ns + nc + nn
        queried <- if (is.data.frame(lookup) && "status" %in% names(lookup)) lookup$status else character()
        gap <- any(queried %in% c("failed", "unavailable", "not-run", "empty"))
        unknown <- length(s) != ncomp + nx + nf
        state <- if (ns && nc) "mixed" else if (ns) "support" else if (nc) "conflict"
            else if (nx && !ncomp) "noncomputable" else "undetermined"
        availability <- if (!valid) "missing-ledger"
            else if (length(s) && (nf || gap || unknown)) {
                if (ncomp || nx) "partial" else if (nf) "failed" else "unknown"
            } else if (length(s)) "available"
            else if (any(queried %in% c("failed", "unavailable"))) "lookup-unavailable"
            else if (any(queried %in% "not-run")) "unqueried"
            else if (length(queried) && all(queried %in% "empty")) "no-dataset"
            else "not-tested"
        data.frame(evidence_status = state, availability = availability,
            incomplete = availability != "available" || nx > 0L, n_attempted = as.integer(length(s)),
            n_computable = as.integer(ncomp), n_support = as.integer(ns), n_conflict = as.integer(nc),
            n_nonsignificant = as.integer(nn), n_noncomputable = as.integer(nx), n_failed = as.integer(nf),
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' @rdname pathway_evidence
#' @export
print.pathway_evidence <- function(x, n = 10L, ...) {
    cat("<pathway evidence>\n")
    evidence <- .pathway_evidence_status(x)
    states <- table(factor(evidence$evidence_status,
        levels = c("support", "conflict", "mixed", "undetermined", "noncomputable")))
    cat("  evidence:     ", paste0(names(states), " ", as.integer(states), collapse = "  "), "\n", sep = "")
    counts <- table(factor(x$level, levels = c("A", "B", "C")))
    cat("  legacy levels (secondary): ", paste0(names(counts), " ", as.integer(counts),
                                  collapse = "  "), "\n", sep = "")
    cat("  Ledger consistency is not a grading pass or causal identification; support can remain when the member-support floor fails.\n")

    ## Preserve legacy grade limits without confusing them with ledger state.
    rs <- x$reason[!is.na(x$reason)]
    if (length(rs)) {
        rc <- sort(table(rs), decreasing = TRUE)
        cat("  legacy capped by: ", paste0(names(rc), " ", as.integer(rc),
                                       collapse = ", "), "\n", sep = "")
    }

    shown <- utils::head(x, n)
    if (nrow(shown)) {
        cat("\n")
        for (i in seq_len(nrow(shown))) {
            row <- shown[i, ]
            driver <- if (!length(row$driver) || is.na(row$driver)) "-" else
                sprintf("%s/%s padj %.2g", row$driver, row$cell_line, row$padj)
            ## `%||%` only catches NULL, not NA; direction used to print as "NA" and support as "NA%".
            dir_txt <- if (is.na(row$direction)) "?" else row$direction
            sup_txt <- if (is.na(row$support)) "-" else sprintf("%.0f%%", 100 * row$support)
            ev <- evidence[i, ]
            cat(sprintf("  [%s] %s (legacy level %s; %s, %d members, %s together)\n",
                        ev$evidence_status, substr(row$pathway, 1, 52), row$level, dir_txt,
                        row$n_members, sup_txt))
            cat(sprintf("      ledger: %d support / %d conflict / %d nonsignificant; %d/%d computable, %d noncomputable, %d failed [%s]\n",
                ev$n_support, ev$n_conflict, ev$n_nonsignificant, ev$n_computable,
                ev$n_attempted, ev$n_noncomputable, ev$n_failed, ev$availability))
            why <- if (length(row$reason) && !is.na(row$reason)) sprintf("  [%s]", row$reason) else ""
            cat(sprintf("      driver: %s%s\n", driver, why))
        }
        if (nrow(x) > nrow(shown)) {
            cat(sprintf("  ... %d more\n", nrow(x) - nrow(shown)))
        }
    }

    ## The statement goes at the end instead of printing on every row: it is a claim to be tested,
    ## not part of the summary.
    a <- x[x$level == "A" & evidence$evidence_status == "support", , drop = FALSE]
    if (nrow(a)) {
        cat("\n  legacy claim to test (not causal identification):\n")
        cat("   ", a$statement[1], "\n")
    }
    invisible(x)
}


## Bring the member IDs in a pathway result into symbol space.
##
## Detection: take the members of every pathway to be reported; if **all** of them are purely
## numeric, treat them as Entrez and translate. Only this one test is made, because outside
## Entrez/RefSeq/Ensembl there is no common purely numeric form, and guessing wrong (silently all
## empty) costs far more than not guessing (asking the caller to convert).
##
## If translation fails, **return unchanged** and let the guard downstream report it -- the error
## message is more complete there.
.pathway_ids_to_symbols <- function(d, signed) {
    if (is.null(d) || !nrow(d) || !"geneID" %in% names(d)) return(d)
    toks <- unique(trimws(unlist(strsplit(as.character(d$geneID), "/"))))
    toks <- toks[nzchar(toks) & !is.na(toks)]
    if (!length(toks) || !all(grepl("^[0-9]+$", toks))) return(d)
    ## Already numeric, and `signed` numeric too: no conversion needed, they are in the same space.
    if (!is.null(signed) && length(signed) && all(grepl("^[0-9]+$", names(signed)))) return(d)
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(d)

    conv <- tryCatch(suppressMessages(
        AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = toks,
                              keytype = "ENTREZID", columns = "SYMBOL")),
        error = function(e) NULL)
    if (is.null(conv) || !nrow(conv)) return(d)
    map <- stats::setNames(as.character(conv$SYMBOL), as.character(conv$ENTREZID))
    d$geneID <- vapply(strsplit(as.character(d$geneID), "/"), function(ids) {
        sym <- unname(map[trimws(ids)])
        paste(stats::na.omit(sym), collapse = "/")
    }, character(1))
    d[nzchar(d$geneID), , drop = FALSE]
}


## When the low grades of a batch of pathways **all** come from one external gap, warn once and say
## so plainly.
##
## Fires only when there is no A at all: if some rows already reached A the conditions are partly
## met, reading `reason` row by row is enough, and one more message is just noise.
.pathway_warn_capped <- function(out, drivers, failed = character(0)) {
    ## `!any(level %in% c("A","B","C"))` is vacuously false: level only takes those three values.
    ## That line could never fire, so it was as good as unwritten.
    if (!nrow(out)) return(invisible(NULL))
    if (any(out$level == "A", na.rm = TRUE)) return(invisible(NULL))

    why <- c("single-cell-line", "lookup-failed", "no-dataset", "too-small", "no-test-completed")
    rs <- out$reason[!is.na(out$reason) & out$reason %in% why]
    if (!length(rs) || length(unique(rs)) != 1L) return(invisible(NULL))
    ## Only report on a majority -- a handful of small pathways cannot explain a whole batch.
    if (length(rs) < max(1L, floor(0.5 * nrow(out)))) return(invisible(NULL))

    ## Name the candidates the sentence is true of. After a partial lookup
    ## failure, a single-cell-line warning that also named the unresolved
    ## candidate claimed it had evidence in one cell line when it had none.
    named <- if (identical(unique(rs), "lookup-failed")) failed else drivers
    if (!length(named)) return(invisible(NULL))
    who <- paste0("`", paste(named, collapse = "`, `"), "`")
    msg <- switch(unique(rs),
        "single-cell-line" = sprintf(
            "level A was not reached: %s has direction-consistent significant evidence in only one cell line for the selected candidate; level A needs two for that same candidate.",
            who),
        "lookup-failed" = sprintf(
            "level A is unreachable: GPSAdb could not be queried for %s. This is a lookup failure, not evidence about the factor.",
            who),
        "no-dataset" = sprintf(
            "level A is unreachable: GPSAdb returned no usable perturbation dataset for %s. The pathways are described, not tested.",
            who),
        "too-small" = paste(
            "level A is unreachable: the pathways were too small for test_gene_set to produce a p-value,",
            "so the perturbation evidence could not be evaluated at all."),
        "no-test-completed" = paste(
            "every GPSAdb call failed, so no perturbation evidence could be collected at all.",
            "The C grades here describe an outage, not the biology -- check the connection and rerun.")
    )
    warning(msg, call. = FALSE)
    invisible(NULL)
}


## Where the pathway name comes from.
##
## The documentation allows two inputs: an `enrichResult` (with `Description` / `ID`) and the
## `$pathways` table from `evidence()` (with `pathway_description` / `pathway`). Recognising only the
## first set of column names gave every row of the second entry point an `NA` -- measured, that is
## exactly what happens, with no error, just an empty pathway-name column in the table.
##
## Both sets of names are tried, and only when neither exists does it fall back to `NA`. This is not
## a guess: `evidence_pathway.R` writes `pathway_description`, which is an **established contract**
## within the same package.
.pathway_label <- function(row, i = 1L) {
    for (nm in c("Description", "pathway_description", "ID", "pathway")) {
        if (!is.null(row[[nm]])) {
            v <- as.character(row[[nm]][1])
            if (!is.na(v) && nzchar(v)) return(v)
        }
    }
    NA_character_
}
