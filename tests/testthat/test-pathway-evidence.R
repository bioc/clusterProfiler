## The first batch of tests for pathway mode.
##
## These test the **decision logic**, not the network. `pathway_evidence()` calls
## GPSAdb when drivers are present, so cases that touch the network are skipped
## outright; the offline cases cover direction inference, member concordance,
## grades and statements.

## Bind several rows into one enrichment result, used for mixed input where one
## entry matches and another does not.
rbind_result <- function(...) {
    parts <- lapply(list(...), function(z) as.data.frame(z, stringsAsFactors = FALSE))
    do.call(rbind, parts)
}

pathway_fixture_result <- function() {
    data.frame(
        ID = c("GO:1", "GO:2", "GO:3"),
        Description = c("response to hypoxia", "coherent pathway", "incoherent pathway"),
        geneID = c("A/B/C/D", "E/F/G/H", "I/J/K/L"),
        Count = 4L, p.adjust = c(1e-10, 1e-5, 1e-3),
        stringsAsFactors = FALSE)
}

test_that("direction is derived from the members when the result has none", {
    ## An ORA result carries no direction of its own: it says the pathway is
    ## enriched, not which way it moved.
    ## The first version returned NULL here, so every row became untestable - that
    ## only showed up once it ran.
    r <- pathway_fixture_result()[1, ]
    members <- c("A", "B", "C", "D")
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5)
    expect_identical(.pathway_direction(r, members, signed), "up")
    expect_identical(.pathway_direction(r, members, -signed), "down")
    ## With no signed, do not guess.
    expect_null(.pathway_direction(r, members, NULL))
})

test_that("member support measures the fraction moving together", {
    ## A high score driven by a few genes is not support for the whole pathway.
    members <- c("A", "B", "C", "D")
    expect_equal(.pathway_support(members, c(A = 1, B = 2, C = 3, D = 4), "up"), 1)
    expect_equal(.pathway_support(members, c(A = 1, B = -2, C = 3, D = -4), "up"), 0.5)
    ## When none of the members is among the tested genes the answer is NA, not 0 -
    ## "not tested" is not "did not pass".
    expect_true(is.na(.pathway_support(members, c(Z = 1), "up")))
})

test_that("no signed input is untestable, and says which", {
    ## No signed means **not tested**, not "tested and refuted".
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, max_pathways = 3)
    expect_equal(nrow(out), 3)
    expect_true(all(out$level == "C"))
    expect_true(all(out$verdict == "untestable"))
    expect_true(all(grepl("No signed fold changes were supplied", out$statement)))
})

test_that("support derived from the members themselves cannot promote to level B", {
    ## A defect found in review: when the direction is inferred from the member
    ## majority, support is vacuously >= 0.5, so `min_support = 0.5` always holds and
    ## every pathway with a direction gets a free B. The direction itself comes from
    ## the members, so using their concordance as evidence has the same genes both
    ## setting the direction and confirming it.
    r <- pathway_fixture_result()[1:2, ]
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = -2, F = -1, G = -3, H = -0.5)
    out <- pathway_evidence(r, signed = signed, max_pathways = 2)
    expect_true(all(out$level == "C"))
    expect_false(any(out$verdict == "coherent-untested"))
    ## But support is still reported as measured: it is not wrong, it just cannot
    ## count as independent evidence.
    expect_true(all(out$support == 1))
})

test_that("an externally supplied direction makes member support meaningful again", {
    ## When the direction comes from the result itself (external information), member
    ## concordance is independent evidence and can earn a B.
    r <- pathway_fixture_result()[1, ]
    r$pathway_direction <- "up"
    signed <- c(A = 2, B = 1.5, C = 3, D = -0.5)
    out <- pathway_evidence(r, signed = signed)
    expect_equal(out$level, "B")
    expect_equal(out$verdict, "coherent-untested")
    expect_equal(out$support, 0.75)
})

test_that("members splitting evenly are unsupported, not untested", {
    ## signed is supplied but the members split evenly: this is **tested, no common
    ## direction**, a different conclusion from "not tested".
    r <- pathway_fixture_result()[3, ]
    signed <- c(I = 2, J = -1, K = 3, L = -4)
    out <- pathway_evidence(r, signed = signed)
    expect_equal(out$level, "C")
    expect_equal(out$verdict, "unsupported")
    expect_match(out$statement, "split evenly")
})

test_that("an empty result returns an empty frame with the right columns", {
    out <- pathway_evidence(pathway_fixture_result()[0, ])
    expect_equal(nrow(out), 0)
    expect_true(all(c("pathway", "level", "statement") %in% names(out)))
})

test_that("level A requires two cell lines, and one dataset is not enough", {
    ## Measured justification: three of four drivers fail when one cell line is
    ## swapped, and the effect-size sign flips too (results/stage3_refutation.tsv).
    ## An implementation that keeps only the strongest dataset makes grade A mean
    ## "one dataset supports it", which is too weak, so grade A now requires
    ## reproduction.
    ##
    ## Only the decision logic can be tested here, with no network call: a skip
    ## guards it and the ground truth is computed inside pathway_evidence.
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))
    r <- pathway_fixture_result()[1, ]
    out <- pathway_evidence(r, drivers = "HIF1A", signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_true("n_cell_lines_supported" %in% names(out))
    expect_true("n_tests_run" %in% names(out))
    ## Significant support from a single cell line cannot be A.
    if (out$n_tests_run > 0 && out$padj < 0.05 && out$n_cell_lines_supported < 2) {
        expect_false(out$level == "A")
        expect_equal(out$verdict, "single-dataset")
    }
})

test_that("every row says why it stopped short of level A", {
    ## Why `reason` exists: `single-dataset` and `direction-conflict` both land in B,
    ## one for too little evidence and one for evidence pointing the other way; the
    ## three C cases (members not found, pathway too small, tested but not
    ## significant) also need separating. Otherwise readers cannot distinguish
    ## "not found" from "cannot be looked for".
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, drivers = NULL, signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_true("reason" %in% names(out))
    ## With no drivers, B can only come from "the pathway moved but no perturbation
    ## test was run".
    if (any(out$level == "B")) {
        expect_true(all(out$reason[out$level == "B"] %in% "no-perturbation-test"))
    }
    ## Grade A should carry no reason for stopping short.
    expect_true(all(is.na(out$reason[out$level == "A"])))
})

test_that("the result carries a class so print.pathway_evidence is used", {
    ## Measured: this code was never dispatched even once after it was written - with
    ## no class the function returns a bare data.frame and the user sees a 20-column
    ## raw table. That is the hardest kind of "looks present, never ran".
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, drivers = NULL, signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_s3_class(out, "pathway_evidence")
    expect_true(is.data.frame(out))
})

test_that("the statement agrees with the numbers printed beside it", {
    ## A defect found in review, and of the **still present after a fix** kind:
    ## `why_c` merged "members move in scattered directions" and "the direction was
    ## inferred from the members themselves" into one `low-support`, so a pathway
    ## with support = 100% printed "100% of members moved up, below the support
    ## floor of 50%". The `else` branch added last round was never reachable.
    ##
    ## The value of this test is that it **asserts the relation between two columns**,
    ## not just their individual values - testing level and support alone would let
    ## the contradiction survive.
    ## It needs >= 5 members, otherwise `too-small` hits first and this branch is
    ## never reached.
    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)
    out <- pathway_evidence(r, signed = signed)

    expect_equal(out$support, 1)
    expect_false(grepl("below the support floor", out$statement))
    expect_match(out$statement, "inferred from these same members")
    expect_equal(out$reason, "member-derived")

    ## The reverse: only a support rate genuinely below the floor may be called
    ## "below the floor".
    r2 <- r
    r2$pathway_direction <- "up"          # external direction, so support is independent evidence
    signed2 <- c(A = 2, B = -1.5, C = -3, D = -0.5, E = -1)
    out2 <- pathway_evidence(r2, signed = signed2)
    expect_lt(out2$support, 0.5)
    expect_match(out2$statement, "below the support floor")
    expect_equal(out2$reason, "low-support")
})

test_that("max_pathways = 0 does not crash", {
    ## There was no guard after `head(d, 0)`: do.call(rbind, list()) returns NULL, and
    ## the next attr(out, ...) <- throws "attempt to set an attribute on NULL".
    r <- pathway_fixture_result()
    out <- expect_no_error(pathway_evidence(r, max_pathways = 0L,
                                            signed = c(A = 2, B = 1.5)))
    ## Clamped to 1 rather than erroring, and not 0 rows either - "give me zero" has
    ## no meaningful answer.
    expect_equal(nrow(out), 1)
    expect_s3_class(out, "pathway_evidence")
})

test_that("an empty result carries the same class as a full one", {
    ## When the two return shapes disagree, an empty result skips the print method and
    ## the user sees two different outputs from the same function.
    expect_s3_class(pathway_evidence(pathway_fixture_result()[0, ]), "pathway_evidence")
})

test_that("print shows a placeholder rather than NA for untestable rows", {
    ## `%||%` catches NULL but not NA, so direction printed as "NA" and support as
    ## "NA%".
    r <- pathway_fixture_result()[1, ]
    out <- pathway_evidence(r)          # no signed -> the direction is NA
    txt <- paste(utils::capture.output(print(out)), collapse = "\n")
    expect_false(grepl("NA%", txt, fixed = TRUE))
    expect_match(txt, "\\?")
})

test_that("a failed perturbation call is reported as an outage, not as weak evidence", {
    skip_if_no_gpsadb()
    ## Found in review: `ask()` treated a call failure and a successful call with no
    ## result both as next, so one outage emptied `tested` and that was then read as
    ## "members disagree on direction". Measured: HIF1A (4 cell lines) under a
    ## simulated outage got `low-support`, with zero warnings.
    ##
    ## The network is needed first to resolve drivers into datasets, so it goes
    ## through the same network switch.
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)

    ns <- asNamespace("gpsadb")
    orig <- get("gpsadb_test_gene_set", envir = ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", orig, ns = "gpsadb"), add = TRUE)
    assignInNamespace("gpsadb_test_gene_set",
                      function(...) stop("simulated outage"), ns = "gpsadb")

    expect_warning(out <- pathway_evidence(r, drivers = "HIF1A", signed = signed),
                   "GPSAdb calls failed")
    expect_equal(out$reason, "no-test-completed")
    expect_equal(out$n_tests_run, 0)
    expect_match(out$statement, "outage")
})

test_that("subsetting columns drops the class so print keeps working", {
    ## `[.data.frame` carries the custom class along, so `out[, c("pathway","level")]`
    ## still claims to be pathway_evidence while print reads driver/reason/support,
    ## columns that no longer exist; `is.na(NULL)` is logical(0) and `if` on length 0
    ## throws "argument is of length zero". Column subsetting is far too common to
    ## leave unhandled.
    r <- pathway_fixture_result()
    out <- pathway_evidence(r, signed = c(A = 2, B = 1.5, C = 3, D = 0.5))

    sub <- out[, c("pathway", "level")]
    expect_false(inherits(sub, "pathway_evidence"))
    expect_no_error(print(sub))

    ## A row subset keeps every column, so the class should stay.
    expect_s3_class(out[1, ], "pathway_evidence")
    ## Subsetting one column yields a vector, which must not be forced to carry the
    ## class.
    expect_false(inherits(out[, "level"], "pathway_evidence"))
})

test_that("empty separators and NA are not counted as members", {
    ## `"A/B/C/D/E/"` used to count 6 members, so `n_members` and the number in the
    ## statement were both one too high; `geneID = NA` produced a member called "NA"
    ## and reported "no members" as one.
    base <- pathway_fixture_result()[1, ]
    mk <- function(g) { base$geneID <- g; base }

    members <- function(g) clusterProfiler:::.pathway_members(mk(g))
    expect_length(members("A/B/C/D/E"), 5)
    expect_length(members(" A/B/C/D/E/ "), 5)
    expect_length(members("A//B/C/D/E"), 5)
    expect_length(clusterProfiler:::.pathway_members(mk(NA_character_)), 0)
})

test_that("no common direction names which of the three cases it is", {
    ## No signed / members absent from the differential table / members genuinely
    ## split evenly - three different conclusions.
    ## The latter two used to share "Members split evenly", so a pathway with 0
    ## members was also called "split evenly".
    base <- pathway_fixture_result()[1, ]

    r1 <- base
    expect_equal(pathway_evidence(r1)$reason, "no-fold-changes")

    r2 <- base
    r2$geneID <- "ZZZ/YYY"
    r2$Description <- "members not in the differential table"
    out2 <- pathway_evidence(rbind_result(r2, base), signed = c(A = 2, B = 1.5, C = 3, D = 0.5))
    expect_equal(out2$reason[out2$pathway == r2$Description], "no-members")
    expect_match(out2$statement[out2$pathway == r2$Description],
                 "appear in the signed gene list")

    r3 <- base
    r3$geneID <- "A/B"
    r3$pathway_direction <- NA
    out3 <- pathway_evidence(r3, signed = c(A = 2, B = -1))
    expect_equal(out3$reason, "no-direction")
    expect_match(out3$statement, "split evenly")
})


test_that("the documented evidence() pathway table keeps its names", {
    ## The documentation allows `x` to be the `$pathways` table of an `evidence()`
    ## bundle, and that table uses the column names `pathway` / `pathway_description`,
    ## not `ID` / `Description`.
    ## Accepting only the latter left every pathway name NA - no error, just empty
    ## names.
    pw <- data.frame(
        pathway = "R-HSA-1234176",
        pathway_description = "Oxygen-dependent proline hydroxylation of HIF Alpha",
        pathway_direction = "up", n_pathway_genes = 5L, pathway_padj = 1e-8,
        regulator = "VHL", source = "reactome", term = "R-HSA-1234176",
        Count = 5L, setSize = 40L, padj_within_pathway = 1e-5,
        geneID = "A/B/C/D/E", stringsAsFactors = FALSE)

    out <- pathway_evidence(pw, signed = c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1))
    expect_false(is.na(out$pathway))
    expect_equal(out$pathway, "Oxygen-dependent proline hydroxylation of HIF Alpha")
})

test_that("a GSEA result takes its direction from NES, which beats the members", {
    ## A `gseaResult` carries NES, so the direction should come from it; the member
    ## majority direction is used only when there is no NES.
    g <- data.frame(ID = "GO:1", Description = "response to hypoxia",
                    setSize = 5L, enrichmentScore = 0.8, NES = -2.0,
                    pvalue = 0.001, p.adjust = 0.01, qvalue = 0.01,
                    rank = 100, leading_edge = "tags=50%",
                    core_enrichment = "A/B/C/D/E", geneID = "A/B/C/D/E",
                    stringsAsFactors = FALSE)
    ## Every member is up-regulated, but NES says the pathway is down - the direction
    ## must follow NES.
    out <- pathway_evidence(g, signed = c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1))
    expect_equal(out$direction, "down")
    ## When the direction opposes the members, support is 0 and cannot count as
    ## evidence.
    expect_equal(out$support, 0)
    expect_equal(out$reason, "low-support")
})

test_that("an outage is not hidden behind 'no perturbation test identified a driver'", {
    skip_if_no_gpsadb()
    ## When the direction is external and support is high enough, the B branch fires
    ## first and reports "every test failed" as "no test could identify a driver".
    ## Measured: an outage produced exactly that sentence, with no way to tell it was a
    ## network problem.
    ##
    ## The network is needed to resolve datasets, so it goes through the same switch.
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)

    ns <- asNamespace("gpsadb")
    orig <- get("gpsadb_test_gene_set", envir = ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", orig, ns = "gpsadb"), add = TRUE)
    assignInNamespace("gpsadb_test_gene_set",
                      function(...) stop("simulated outage"), ns = "gpsadb")

    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A", signed = signed))
    expect_equal(out$reason, "no-test-completed")
    expect_equal(out$level, "C")
})

test_that("a driver that was tested and refuted is not called untested", {
    skip_if_no_gpsadb()
    ## Tests ran and were all non-significant, versus never ran at all: two different
    ## things. Both used to share `coherent-untested`, so a pathway tested in 4
    ## datasets with no significance was called "untested" - the word `untested`
    ## itself was wrong.
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("CLUSTERPROFILER_NETWORK_TESTS"), "1"))

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    signed <- c(A = 2, B = 1.5, C = 3, D = 0.5, E = 1)

    ns <- asNamespace("gpsadb")
    orig <- get("gpsadb_test_gene_set", envir = ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", orig, ns = "gpsadb"), add = TRUE)
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 0.5, es = -0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")

    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A", signed = signed))
    expect_gt(out$n_tests_run, 0)
    expect_equal(out$verdict, "coherent-refuted")
    expect_equal(out$reason, "not-significant")
    expect_match(out$statement, "found no significant movement")

    ## With no drivers it is still "not tested"; the two must not be conflated.
    out2 <- pathway_evidence(r, signed = signed)
    expect_equal(out2$verdict, "coherent-untested")
    expect_equal(out2$reason, "no-perturbation-test")
})

test_that("the level-A statement admits what it does not establish", {
    skip_if_no_gpsadb()
    ## Measured: on the same hypoxia pathway, taking 30 drivers at random from the
    ## inventory, 4 of them (all essential genes) also got an A. So A only says
    ## "knocking it down moves this pathway", not "it is the driver".
    ## The statement must contain that sentence, otherwise readers take "consistent"
    ## for "identified".
    ##
    ## Stubbing out a genuine grade A (two cell lines, same direction, significant),
    ## without the network.
    ## The earlier version used a 3-member pathway that could never reach A, so the
    ## whole test was skipped - skipping is not coverage, hence the stubs.
    ns <- asNamespace("gpsadb")
    orig_s <- get("gpsadb_search_datasets", envir = ns)
    orig_t <- get("gpsadb_test_gene_set", envir = ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", orig_s, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", orig_t, ns = "gpsadb")
    }, add = TRUE)

    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        list(data.frame(perturbed_gene = rep(perturbed_gene, 2),
                        method = "knockdown",
                        cell_line = c("HeLa", "A-498"),
                        n_samples = c(4L, 4L), in_reference = c(TRUE, TRUE),
                        dataset_id = c("D001", "D002"), stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 1e-6, es = -0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"          # direction from outside; es < 0 means activator
    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A",
                                             signed = c(A = 2, B = 1.5, C = 3,
                                                        D = 0.5, E = 1)))

    expect_equal(out$level, "A")
    expect_match(out$statement, "Consistent with drive")
    expect_match(out$statement, "not evidence that")
})

test_that("a single-cell-line row says the cap is the data, not the evidence", {
    skip_if_no_gpsadb()
    ## There are two kinds of B: direction conflict (genuinely opposed) and single
    ## cell line (merely not enough).
    ## Measured, the latter pointed correctly on all 11 datasets with known ground
    ## truth, so the wording must not describe it as defective.
    ns <- asNamespace("gpsadb")
    orig_s <- get("gpsadb_search_datasets", envir = ns)
    orig_t <- get("gpsadb_test_gene_set", envir = ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", orig_s, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", orig_t, ns = "gpsadb")
    }, add = TRUE)
    ## Only one cell line - exactly the kind of data that cannot reach A.
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        list(data.frame(perturbed_gene = perturbed_gene, method = "knockdown",
                        cell_line = "HeLa", n_samples = 4L, in_reference = TRUE,
                        dataset_id = "D001", stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 1e-6, es = -0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")

    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A",
                                             signed = c(A = 2, B = 1.5, C = 3, D = 2, E = 1)))
    expect_equal(out$level[1], "B")
    expect_equal(out$reason[1], "single-cell-line")
    expect_match(out$statement[1], "limit of the available data")
    expect_match(out$statement[1], "consistent perturbation evidence, but under-evidenced")
})

test_that("a direction conflict is not described as thin data", {
    skip_if_no_gpsadb()
    ## The two kinds of B must be worded apart: a direction conflict is evidence
    ## against, a single cell line is not enough data.
    ns <- asNamespace("gpsadb")
    orig_s <- get("gpsadb_search_datasets", envir = ns)
    orig_t <- get("gpsadb_test_gene_set", envir = ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", orig_s, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", orig_t, ns = "gpsadb")
    }, add = TRUE)
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        list(data.frame(perturbed_gene = perturbed_gene, method = "knockdown",
                        cell_line = "HeLa", n_samples = 4L, in_reference = TRUE,
                        dataset_id = "D001", stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    ## The pathway is up and es is > 0 too - knocking it down pushes it higher, which
    ## reads as repression rather than driving.
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        list(data.frame(dataset_id = id, n_genes = length(genes),
                        n_in_dataset = length(genes), padj = 1e-6, es = 0.9,
                        stringsAsFactors = FALSE))
    }, ns = "gpsadb")
    r <- pathway_fixture_result()[1, ]
    r$geneID <- "A/B/C/D/E"
    r$pathway_direction <- "up"
    out <- suppressWarnings(pathway_evidence(r, drivers = "HIF1A",
                                             signed = c(A = 2, B = 1.5, C = 3, D = 2, E = 1)))
    expect_equal(out$verdict[1], "direction-conflict")
    expect_match(out$statement[1], "repression")
    expect_false(grepl("limit of the available data", out$statement[1]))
})

## Namespace stubs are restored even when an assertion fails; no network calls.
with_pathway_stubs <- function(datasets, outcomes, code) {
    skip_if_no_gpsadb()
    ns <- asNamespace("gpsadb")
    old_search <- get("gpsadb_search_datasets", ns)
    old_test <- get("gpsadb_test_gene_set", ns)
    on.exit({
        assignInNamespace("gpsadb_search_datasets", old_search, ns = "gpsadb")
        assignInNamespace("gpsadb_test_gene_set", old_test, ns = "gpsadb")
    }, add = TRUE)
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        z <- datasets[datasets$driver == perturbed_gene, , drop = FALSE]
        z$n_samples <- rep(4L, nrow(z))
        z$in_reference <- rep(TRUE, nrow(z))
        list(z)
    }, ns = "gpsadb")
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) {
        z <- outcomes[[id]]
        if (inherits(z, "error")) stop(z)
        list(z)
    }, ns = "gpsadb")
    force(code)
}

pathway_test_input <- function() {
    data.frame(Description = "programme", geneID = "A/B/C/D/E",
               pathway_direction = "up", stringsAsFactors = FALSE)
}

pathway_test_signed <- c(A = 1, B = 2, C = 3, D = 4, E = 5)

test_that("replication never pools cell lines across candidates", {
    ds <- data.frame(driver = c("TF1", "TF2"), dataset_id = c("D1", "D2"),
                     cell_line = c("cellA", "cellB"))
    got <- list(D1 = data.frame(es = -1, padj = 0.001),
                D2 = data.frame(es = -1, padj = 0.002))
    with_pathway_stubs(ds, got, {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(),
            drivers = c("TF1", "TF2"), signed = pathway_test_signed))
        expect_equal(out$level, "B")
        expect_equal(out$n_cell_lines_supported, 1L)
        expect_equal(out$driver, "TF1")
        ledger <- attr(out, "perturbation_tests")
        expect_equal(ledger$driver, c("TF1", "TF2"))
        expect_equal(ledger$status, c("supports", "supports"))
    })
})

test_that("significance cannot bypass measurable member support", {
    ds <- data.frame(driver = "TF1", dataset_id = c("D1", "D2"),
                     cell_line = c("cellA", "cellB"))
    got <- list(D1 = data.frame(es = -1, padj = 0.001),
                D2 = data.frame(es = -1, padj = 0.002))
    with_pathway_stubs(ds, got, {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = c(A = 1, B = -2, C = -3, D = -4, E = -5)))
        expect_equal(out$level, "C")
        expect_equal(out$reason, "low-support")
        expect_match(out$statement, "below the support floor")
        expect_equal(attr(out, "perturbation_tests")$status, rep("supports", 2))
        inferred <- pathway_test_input()
        inferred$pathway_direction <- NULL
        strict <- suppressWarnings(pathway_evidence(inferred, drivers = "TF1",
            min_support = 0.8, signed = c(A = 1, B = 2, C = 3, D = -4, E = -5)))
        expect_equal(strict$level, "C")
        expect_equal(strict$reason, "low-support")
        unknown <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1"))
        expect_equal(unknown$level, "C")
        ## Direction is given, but no fold changes: the missing input is the
        ## measurements, not the identifiers.
        expect_equal(unknown$reason, "no-fold-changes")
    })
})

test_that("the ledger retains all outcomes and warns about partial failures", {
    ds <- data.frame(driver = "TF1", dataset_id = paste0("D", 1:6),
                     cell_line = paste0("cell", 1:6))
    got <- list(D1 = data.frame(es = -1, padj = 0.001),
                D2 = data.frame(es = -1, padj = 0.002),
                D3 = data.frame(es = 1, padj = 0.003),
                D4 = data.frame(es = -1, padj = 0.5),
                D5 = simpleError("offline failure"),
                D6 = data.frame(es = NA_real_, padj = NA_real_))
    with_pathway_stubs(ds, got, {
        expect_warning(out <- pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = pathway_test_signed), "1 of 6 GPSAdb calls failed")
        ledger <- attr(out, "perturbation_tests")
        expect_identical(ledger$status, c("supports", "supports", "conflicts",
                                         "not-significant", "failed", "not-computable"))
        expect_identical(ledger$direction, rep("up", 6))
        expect_type(ledger$es, "double")
        expect_type(ledger$padj, "double")
        expect_equal(out$n_tests_run, 4L)
        expect_equal(out$n_cell_lines_supported, 2L)
        expect_equal(out$level, "A")
        expect_match(out$statement, "contradictory across datasets")
        expect_equal(attr(out, "network_calls"), 6L)
    })
})

test_that("empty and untestable results retain a typed ledger", {
    empty <- suppressWarnings(pathway_evidence(pathway_test_input()[0, ]))
    ledger <- attr(empty, "perturbation_tests")
    expect_identical(names(ledger), c("pathway", "driver", "dataset", "cell_line",
                                     "es", "padj", "status", "direction"))
    expect_equal(nrow(ledger), 0L)
    expect_type(ledger$es, "double")
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = simpleError("must not call")), {
        x <- pathway_test_input()
        x$pathway_direction <- NULL
        out <- suppressWarnings(pathway_evidence(x, drivers = "TF1"))
        expect_equal(attr(out, "network_calls"), 0L)
        expect_equal(attr(out, "perturbation_tests")$status, "not-computable")
    })
})

test_that("memoisation preserves each pathway outcome and outage status", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = simpleError("offline failure")), {
        x <- rbind(pathway_test_input(), pathway_test_input())
        x$Description <- c("first", "second")
        warnings <- character()
        out <- withCallingHandlers(pathway_evidence(x, drivers = "TF1",
            signed = pathway_test_signed), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        expect_true(any(grepl("All 1 GPSAdb calls failed", warnings)))
        expect_equal(attr(out, "network_calls"), 1L)
        expect_equal(attr(out, "perturbation_tests")$pathway, c("first", "second"))
        expect_equal(attr(out, "perturbation_tests")$status, rep("failed", 2))
        expect_equal(out$reason, rep("no-test-completed", 2))
    })
})

test_that("non-significance is not described as proof of falsity", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = data.frame(es = -1, padj = 0.5)), {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = pathway_test_signed))
        expect_equal(out$verdict, "coherent-refuted") # legacy API key
        expect_equal(out$reason, "not-significant")
        expect_equal(attr(out, "perturbation_tests")$status, "not-significant")
        expect_match(out$statement, "found no significant movement")
    })
})

test_that("requested candidates survive unresolved datasets and empty results", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    requested <- c("TF1", "UNRESOLVED")
    with_pathway_stubs(ds, list(D1 = data.frame(es = -1, padj = 0.001)), {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(),
            drivers = requested, signed = pathway_test_signed))
        expect_identical(attr(out, "drivers"), requested)
        empty <- suppressWarnings(pathway_evidence(pathway_test_input()[0, ], drivers = requested))
        expect_identical(attr(empty, "drivers"), requested)
        expect_null(attr(.pathway_evidence_empty(), "drivers"))
    })
})

test_that("row subsets scope the ledger without losing requested candidates", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    with_pathway_stubs(ds, list(D1 = data.frame(es = -1, padj = 0.001)), {
        x <- rbind(pathway_test_input(), pathway_test_input())
        x$Description <- c("first", "second")
        out <- suppressWarnings(pathway_evidence(x, drivers = c("TF1", "UNRESOLVED"),
            signed = pathway_test_signed))
        expect_equal(attr(out[1, ], "perturbation_tests")$pathway, "first")
        expect_equal(attr(out[2, , drop = FALSE], "perturbation_tests")$pathway, "second")
        expect_equal(nrow(attr(out[FALSE, ], "perturbation_tests")), 0L)
        expect_identical(attr(out[FALSE, ], "drivers"), c("TF1", "UNRESOLVED"))
        expect_equal(attr(out[1, "level", drop = FALSE], "perturbation_tests")$pathway, "first")
        expect_equal(nrow(attr(out["pathway"], "perturbation_tests")), 2L)
    })
})

test_that("blank missing or padded cell lines cannot fabricate replication", {
    for (other in c("", " ", NA_character_, " cellA ")) {
        ds <- data.frame(driver = "TF1", dataset_id = c("D1", "D2"),
                         cell_line = c("cellA", other))
        got <- list(D1 = data.frame(es = -1, padj = 0.001),
                    D2 = data.frame(es = -1, padj = 0.002))
        with_pathway_stubs(ds, got, {
            warnings <- character()
            out <- withCallingHandlers(pathway_evidence(pathway_test_input(), drivers = "TF1",
                signed = pathway_test_signed), warning = function(w) {
                    warnings <<- c(warnings, conditionMessage(w))
                    invokeRestart("muffleWarning")
                })
            expect_equal(out$level, "B")
            expect_equal(out$n_cell_lines_supported, 1L)
            expect_true(any(grepl("Level A is unreachable", warnings)))
        })
    }
})

test_that("conflicting identities for a dataset key are rejected", {
    ds <- data.frame(driver = "TF1", dataset_id = c("D1", "D1", "D2"),
                     cell_line = c("cellA", "cellB", "cellC"))
    got <- list(D1 = simpleError("ambiguous dataset must not run"),
                D2 = data.frame(es = -1, padj = 0.001))
    with_pathway_stubs(ds, got, {
        warnings <- character()
        out <- withCallingHandlers(pathway_evidence(pathway_test_input(), drivers = "TF1",
            signed = pathway_test_signed), warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        expect_true(any(grepl("conflicting cell-line identities: D1", warnings)))
        expect_equal(out$level, "B")
        expect_equal(out$n_cell_lines_supported, 1L)
        expect_equal(attr(out, "network_calls"), 1L)
        expect_equal(attr(out, "perturbation_tests")$dataset, "D2")
    })
})

test_that("Entrez conversion includes IDs beyond the fifth pathway", {
    skip_if_not_installed("org.Hs.eg.db")
    ns <- asNamespace("AnnotationDbi")
    original <- get("select", ns)
    on.exit(assignInNamespace("select", original, ns = "AnnotationDbi"), add = TRUE)
    seen <- NULL
    assignInNamespace("select", function(x, keys, ...) {
        seen <<- keys
        data.frame(ENTREZID = keys, SYMBOL = paste0("G", keys))
    }, ns = "AnnotationDbi")
    x <- data.frame(Description = paste0("path", 1:6), geneID = as.character(1:6))
    out <- .pathway_ids_to_symbols(x, c(G1 = 1))
    expect_identical(seen, as.character(1:6))
    expect_equal(nrow(out), 6L)
    expect_identical(out$geneID, paste0("G", 1:6))
})

test_that("dataset resolution accepts every removal method, not only 'knockdown'", {
    ## GPSAdb spells removal siRNA/shRNA/CRISPRi/knockout. A literal
    ## method = "knockdown" filter reported STAT1 -- CRISPRi plus knockout,
    ## no literal knockdown -- as having no data, the same false negative the
    ## inventory lookup was fixed for. Overexpression stays excluded.
    d <- data.frame(
        dataset_id = c("d1", "d2", "d3", "d4", "d5"),
        cell_line = c("A", "B", "C", "D", "E"),
        method = c("CRISPRi", "knockout", "siRNA", "overexpression", "shRNA"),
        in_reference = c(TRUE, TRUE, TRUE, TRUE, FALSE),
        n_samples = c(3, 2, 4, 9, 1), stringsAsFactors = FALSE)
    got <- .pathway_usable_datasets(d)
    expect_identical(got$dataset_id, c("d1", "d2", "d3"))
    expect_false(any(grepl("overexpress", got$method)))
})

test_that("unusable dataset tables are empty, not fabricated data", {
    expect_null(.pathway_usable_datasets(NULL))
    expect_null(.pathway_usable_datasets(data.frame()))
    expect_null(.pathway_usable_datasets(data.frame(dataset_id = "d1")))
    d <- data.frame(dataset_id = c("d1", "d1"), cell_line = c("A", " "),
                    method = "siRNA", in_reference = TRUE, stringsAsFactors = FALSE)
    expect_null(suppressWarnings(.pathway_usable_datasets(d)))
    blank <- data.frame(dataset_id = c("d1", "d2"), cell_line = c("A", "  "),
                        method = "siRNA", in_reference = TRUE, stringsAsFactors = FALSE)
    out <- .pathway_usable_datasets(blank)
    expect_equal(nrow(out), 2L)
    expect_true(is.na(out$cell_line[2]))
})

test_that("a failed dataset lookup is not reported as absent data", {
    skip_if_no_gpsadb()
    ns <- asNamespace("gpsadb")
    original <- get("gpsadb_search_datasets", ns)
    on.exit(assignInNamespace("gpsadb_search_datasets", original, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_search_datasets",
                      function(...) stop("simulated outage"), ns = "gpsadb")
    warnings <- character()
    out <- withCallingHandlers(
        pathway_evidence(pathway_test_input(), drivers = "TF1",
                         signed = pathway_test_signed),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_equal(out$reason, "lookup-failed")
    expect_match(out$statement, "lookup failure, not evidence")
    expect_true(any(grepl("lookup failure", warnings)))
    expect_false(any(grepl("No usable perturbation data found", warnings)))
    lookup <- attr(out, "dataset_lookup")
    expect_equal(lookup$driver, "TF1")
    expect_equal(lookup$status, "failed")
})

test_that("an empty lookup stays distinguishable from a failed one", {
    skip_if_no_gpsadb()
    ns <- asNamespace("gpsadb")
    original <- get("gpsadb_search_datasets", ns)
    on.exit(assignInNamespace("gpsadb_search_datasets", original, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_search_datasets",
                      function(...) list(data.frame(dataset_id = character(),
                          cell_line = character(), method = character(),
                          in_reference = logical(), n_samples = integer())),
                      ns = "gpsadb")
    warnings <- character()
    out <- withCallingHandlers(
        pathway_evidence(pathway_test_input(), drivers = "TF1",
                         signed = pathway_test_signed),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_equal(out$reason, "no-dataset")
    expect_true(any(grepl("No usable perturbation data found", warnings)))
    expect_false(any(grepl("lookup failure", warnings)))
    expect_equal(attr(out, "dataset_lookup")$status, "empty")
})

test_that("a partial lookup failure is disclosed even when another driver resolves", {
    skip_if_no_gpsadb()
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    got <- list(D1 = data.frame(es = -1, padj = 0.001))
    ns <- asNamespace("gpsadb")
    original <- get("gpsadb_search_datasets", ns)
    on.exit(assignInNamespace("gpsadb_search_datasets", original, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_search_datasets", function(perturbed_gene, ...) {
        if (identical(perturbed_gene, "TF2")) stop("simulated outage for TF2")
        z <- ds[ds$driver == perturbed_gene, , drop = FALSE]
        z$n_samples <- 4L
        z$in_reference <- TRUE
        list(z)
    }, ns = "gpsadb")
    original_test <- get("gpsadb_test_gene_set", ns)
    on.exit(assignInNamespace("gpsadb_test_gene_set", original_test, ns = "gpsadb"),
            add = TRUE)
    assignInNamespace("gpsadb_test_gene_set", function(id, genes) list(got[[id]]),
                      ns = "gpsadb")
    warnings <- character()
    out <- withCallingHandlers(
        pathway_evidence(pathway_test_input(), drivers = c("TF1", "TF2"),
                         signed = pathway_test_signed),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_true(any(grepl("could not be queried for `TF2`", warnings)))
    expect_true(any(grepl("lookup failure", warnings)))
    ## The single-cell-line warning must name the resolved candidate only: TF2
    ## was never queried, so it cannot have evidence in one cell line.
    single <- warnings[grepl("single cell line", warnings)]
    expect_true(length(single) >= 1L)
    expect_true(any(grepl("TF1", single)))
    expect_false(any(grepl("TF2", single)))
    lookup <- attr(out, "dataset_lookup")
    expect_equal(lookup$status[lookup$driver == "TF1"], "available")
    expect_equal(lookup$status[lookup$driver == "TF2"], "failed")
})

test_that("missing fold changes are not blamed on the identifier space", {
    ## A documented call: direction from pathway_direction/NES, no `signed`.
    ## The perturbation tests run and support the direction, but member support
    ## is unknown, so the row cannot pass the support floor. It must say the
    ## fold changes are missing, not that the members are.
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    got <- list(D1 = data.frame(es = -1, padj = 0.001))
    with_pathway_stubs(ds, got, {
        out <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1"))
        expect_equal(out$level, "C")
        expect_equal(out$reason, "no-fold-changes")
        expect_match(out$statement, "No signed fold changes")
        expect_false(grepl("signed gene list", out$statement))
        expect_true(is.na(out$support))
        expect_equal(attr(out, "perturbation_tests")$status, "supports")
    })
})

test_that("pathway presentation uses ledger evidence even below the support floor", {
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    got <- list(D1 = data.frame(es = -1, padj = 0.001))
    with_pathway_stubs(ds, got, {
        signed <- pathway_test_signed
        signed[1:2] <- -abs(signed[1:2])
        pe <- suppressWarnings(pathway_evidence(pathway_test_input(), drivers = "TF1",
            min_support = 0.8, signed = signed))
        expect_equal(pe$reason, "low-support")
        expect_equal(clusterProfiler:::.pathway_evidence_status(pe)$evidence_status, "support")
        printed <- paste(capture.output(print(pe)), collapse = "\n")
        expect_match(printed, "\\[support\\].*legacy level C")
        expect_match(printed, "not a grading pass")
        expect_match(printed, "1/1 computable")
        attr(pe, "perturbation_tests") <- NULL
        pe$level <- "A"
        printed <- paste(capture.output(print(pe)), collapse = "\n")
        expect_match(printed, "\\[undetermined\\].*legacy level A")
        expect_match(printed, "missing-ledger")
        expect_false(grepl("claim to test:", printed, fixed = TRUE))
    })
    empty <- clusterProfiler:::.pathway_evidence_empty()
    states <- clusterProfiler:::.pathway_evidence_status(empty)
    expect_equal(nrow(states), 0L)
    expect_type(states$evidence_status, "character")
    expect_type(states$n_support, "integer")
    expect_no_error(print(empty))
})

test_that("a column-only subset keeps the retained ledger", {
    two <- data.frame(Description = c("p1", "p2"), geneID = "A/B/C/D/E",
                      pathway_direction = "up", stringsAsFactors = FALSE)
    ds <- data.frame(driver = "TF1", dataset_id = "D1", cell_line = "cellA")
    got <- list(D1 = data.frame(es = -1, padj = 0.001))
    with_pathway_stubs(ds, got, {
        pe <- suppressWarnings(pathway_evidence(two, drivers = "TF1",
            signed = pathway_test_signed))
        reordered <- pe[, rev(names(pe))]
        expect_s3_class(reordered, "pathway_evidence")
        expect_equal(nrow(attr(reordered, "perturbation_tests")), 2L)
        expect_identical(attr(reordered, "perturbation_tests"),
                         attr(pe, "perturbation_tests"))
        rep <- mechanism_report(reordered, inventory = data.frame(perturbed_gene = "TF1"))
        expect_equal(nrow(rep$perturbation_tests), 2L)
        expect_equal(nrow(rep$candidate_evidence), 2L)
        expect_equal(rep$verdicts$evidence_status, c("support", "support"))
        expect_match(paste(capture.output(print(reordered)), collapse = "\n"), "\\[support\\] p1 \\(legacy level B")
        ## Dropping a required column still falls back to a plain data.frame.
        slim <- pe[, c("pathway", "level")]
        expect_false(inherits(slim, "pathway_evidence"))
        slim_report <- mechanism_report(slim, inventory = data.frame(perturbed_gene = "TF1"))
        expect_equal(slim_report$verdicts$evidence_status, c("support", "support"))
        expect_equal(clusterProfiler:::.pathway_evidence_status(pe[2, ])$n_support, 1L)
    })
})
