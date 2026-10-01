library(testthat)
library(clusterProfiler)

## Helper: check that a function's formals include eps with the expected default
check_eps_in_formals <- function(fn, default_val = 1e-10) {
    fmls <- formals(fn)
    expect_true("eps" %in% names(fmls),
                info = sprintf("%s should have 'eps' in its formals", deparse(substitute(fn))))
    expect_equal(fmls$eps, default_val,
                 info = sprintf("%s$eps should default to %s", deparse(substitute(fn)), default_val))
}

## Helper: check that ... is present
check_dots_in_formals <- function(fn) {
    fmls <- formals(fn)
    expect_true("..." %in% names(fmls),
                info = sprintf("%s should have '...' in its formals", deparse(substitute(fn))))
}

# ---- Signature tests ----

test_that("GSEA has eps parameter with correct default", {
    check_eps_in_formals(GSEA)
    check_dots_in_formals(GSEA)
})

test_that("gseGO has eps parameter with correct default", {
    check_eps_in_formals(gseGO)
    check_dots_in_formals(gseGO)
})

test_that("gseMKEGG has eps parameter with correct default", {
    check_eps_in_formals(gseMKEGG)
    check_dots_in_formals(gseMKEGG)
})

test_that("gseKEGG has eps parameter with correct default", {
    check_eps_in_formals(gseKEGG)
    check_dots_in_formals(gseKEGG)
})

# ---- Parameter forwarding tests via body inspection ----

## Verify that the function body of each wrapper actually passes eps and ...
## to enrichit::gsea_gson()
check_body_forwards_eps <- function(fn) {
    fn_name <- deparse(substitute(fn))
    body_text <- deparse(body(fn))
    combined <- paste(body_text, collapse = "\n")

    expect_match(combined, "eps\\s*=\\s*eps",
                 info = sprintf("%s body should forward eps = eps", fn_name))
    expect_match(combined, "\\.\\.\\.",
                 info = sprintf("%s body should forward ...", fn_name))
}

test_that("GSEA body forwards eps and ... to enrichit::gsea_gson", {
    check_body_forwards_eps(GSEA)
})

test_that("gseGO body forwards eps and ... to enrichit::gsea_gson", {
    check_body_forwards_eps(gseGO)
})

test_that("gseMKEGG body forwards eps and ... to enrichit::gsea_gson", {
    check_body_forwards_eps(gseMKEGG)
})

test_that("gseKEGG body forwards eps and ... to enrichit::gsea_gson", {
    check_body_forwards_eps(gseKEGG)
})

# ---- Behavioral test: scoreType forwarding via gseGO ----
## This is the core issue #44 test: when all geneList values are positive
## and scoreType = "pos", there should be no warning about scoreType being "std".

test_that("gseGO forwards scoreType='pos' without spurious warning (issue #44)", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("enrichit")

    library(org.Hs.eg.db)

    set.seed(42)
    all_positive_genes <- sort(
        setNames(
            runif(100, 0.01, 3),
            sample(keys(org.Hs.eg.db, keytype = "SYMBOL"), 100)
        ),
        decreasing = TRUE
    )

    ## Capture warnings
    w <- NULL
    result <- withCallingHandlers(
        gseGO(
            geneList    = all_positive_genes,
            OrgDb       = org.Hs.eg.db,
            keyType     = "SYMBOL",
            ont         = "BP",
            minGSSize   = 5,
            maxGSSize   = 500,
            scoreType   = "pos",
            seed        = 1,
            verbose     = FALSE
        ),
        warning = function(cond) {
            w <<- c(w, list(conditionMessage(cond)))
            invokeRestart("muffleWarning")
        }
    )

    ## The key assertion: no warning mentioning scoreType being "std"
    scoreType_warnings <- Filter(function(msg) {
        grepl("scoreType", msg, ignore.case = TRUE) &&
        grepl("std", msg, ignore.case = TRUE)
    }, w)

    expect_equal(length(scoreType_warnings), 0,
                 info = "gseGO with scoreType='pos' should not warn about scoreType being 'std'")
})

# ---- Behavioral test: seed forwarding via gseGO ----

test_that("gseGO forwards seed parameter for reproducible results", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("enrichit")

    library(org.Hs.eg.db)

    set.seed(42)
    genes <- sort(
        setNames(
            runif(100, 0.01, 3),
            sample(keys(org.Hs.eg.db, keytype = "SYMBOL"), 100)
        ),
        decreasing = TRUE
    )

    ## Run twice with same seed and scoreType - should produce identical results
    res1 <- gseGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "BP", minGSSize = 5, maxGSSize = 500,
                   scoreType = "pos", seed = 123, verbose = FALSE)
    res2 <- gseGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "BP", minGSSize = 5, maxGSSize = 500,
                   scoreType = "pos", seed = 123, verbose = FALSE)

    expect_identical(as.data.frame(res1), as.data.frame(res2))
})

# ---- Behavioral test: eps forwarding via GSEA with TERM2GENE ----

test_that("GSEA forwards eps parameter without error", {
    skip_if_not_installed("enrichit")

    set.seed(42)
    geneList <- sort(setNames(runif(200, 0.01, 3), paste0("gene", 1:200)),
                     decreasing = TRUE)

    TERM2GENE <- data.frame(
        term = rep(paste0("PATH", 1:10), each = 20),
        gene = paste0("gene", 1:200),
        stringsAsFactors = FALSE
    )

    ## Should not error with explicit eps
    expect_no_error(
        GSEA(geneList, minGSSize = 5, maxGSSize = 50,
             eps = 1e-10, scoreType = "pos", TERM2GENE = TERM2GENE, verbose = FALSE)
    )
})
