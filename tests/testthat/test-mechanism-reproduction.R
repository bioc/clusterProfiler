## Reproduction rate (mechanism_reproduction) and its use inside pathway_evidence.
##
## The tests at this layer care only about the **decision rules** and touch no
## network: the input is a constructed matrix and the output a predictable grade.
## The reproduction rate of real data is computed by the offline scripts in `scripts/`.

test_that("concordance is measured against the majority, not against the pathway", {
    ## An activating driver that is knocked down -> members should be
    ## **down-regulated**. Judging by "agrees with the pathway direction" gives a
    ## true driver a reproduction of 0 (HIF1A in practice: 0.00 under the old rule,
    ## 1.00 under the new one).
    profiles <- matrix(0, nrow = 6, ncol = 3,
                       dimnames = list(paste0("g", 1:6), c("d1", "d2", "d3")))
    profiles[1:4, ] <- -1        # three quarters of the members go down, all three datasets agree
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)))
    expect_equal(out$concordance, 1)
    expect_equal(out$majority_sign, -1L)
    expect_equal(out$grade, "reproduces")
})

test_that("a single dataset is not-computable rather than perfectly concordant", {
    ## At k = 1 concordance is vacuously 1. Without a special case, "only one
    ## dataset" is scored as "perfect reproduction", which turns the case that most
    ## needs caution into the strongest conclusion.
    profiles <- matrix(-1, nrow = 4, ncol = 1,
                       dimnames = list(paste0("g", 1:4), "d1"))
    meta <- data.frame(dataset = "d1", driver = "TF1", biosample = "A")
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)))
    expect_true(is.na(out$concordance))
    expect_equal(out$grade, "not-computable")
})

test_that("runs that disagree are graded not-reproduced", {
    profiles <- matrix(c(-1, -1, 1), nrow = 1, ncol = 3,
                       dimnames = list("g1", c("d1", "d2", "d3")))
    profiles <- rbind(profiles, g2 = c(-1, 1, 1), g3 = c(-1, 1, 1))
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(profiles, meta, list(pw = c("g1", "g2", "g3")))
    ## Member means: d1 = -1, d2 = +1/3, d3 = +1/3 -> the majority is positive,
    ## concordance 2/3
    expect_equal(out$concordance, 2 / 3)
    expect_equal(out$grade, "consistent")
})

test_that("max_log2fc drops weak knockdowns before counting datasets", {
    ## Measured: every ARNT / ESR1 / NR3C1 experiment in KnockTF is weaker than -1;
    ## unfiltered they get 4 / 6 / 1 datasets, filtered they drop to zero. The caller
    ## must be able to see this difference.
    profiles <- matrix(-1, nrow = 4, ncol = 3,
                       dimnames = list(paste0("g", 1:4), c("d1", "d2", "d3")))
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"),
                       logFC = c(-2, -0.2, -0.3))
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)),
                                  max_log2fc = -1)
    expect_equal(out$k, 1L)
    expect_equal(out$grade, "not-computable")
})

test_that("reproduces needs three datasets from at least two biosamples", {
    profiles <- matrix(-1, nrow = 4, ncol = 3,
                       dimnames = list(paste0("g", 1:4), c("d1", "d2", "d3")))
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "A", "A"))
    out <- mechanism_reproduction(profiles, meta, list(pw = paste0("g", 1:4)))
    expect_equal(out$concordance, 1)
    ## Three datasets from one cell line are not independent reproduction
    expect_equal(out$grade, "consistent")
})

test_that("the level-B wording has three distinguishable branches", {
    ## No network: the wording helper is tested directly. The three branches must be
    ## distinguishable - if "no second dataset" and "not reproduced" shared one
    ## sentence, readers could not tell "nobody checked" from "checked and failed".
    base <- data.frame(driver = "TF1", pathway = "pw", k = 2L, n_biosamples = 2L,
                       n_genes = 3L, concordance = 1, majority_sign = -1L,
                       grade = "reproduces", stringsAsFactors = FALSE)
    unrep <- base; unrep$concordance <- 0; unrep$grade <- "not-reproduced"
    nocov <- base; nocov$k <- 1L; nocov$concordance <- NA_real_
    nocov$grade <- "not-computable"

    st_rep <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, base)
    expect_match(st_rep, "reproduces across 2 independent datasets", fixed = TRUE)
    expect_match(st_rep, "2/2 agree", fixed = TRUE)

    st_unrep <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, unrep)
    expect_match(st_unrep, "unverified", fixed = TRUE)
    expect_match(st_unrep, "0/2 agree", fixed = TRUE)

    st_nocov <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, nocov)
    expect_match(st_nocov, "No independent reproduction is available", fixed = TRUE)

    ## Without a reproduction table, fall back to the original wording
    st_none <- .pe_single_line_statement("TF1", 3L, "D1", 0.01, NULL)
    expect_match(st_none, "limit of the available data", fixed = TRUE)

    expect_false(identical(st_rep, st_unrep))
    expect_false(identical(st_unrep, st_nocov))
})

test_that("the majority sign is the mode of signs, not the sign of the sum", {
    ## One dataset with a very large effect size must not overturn two smaller but
    ## agreeing ones. The old implementation used sign(sum(shift)):
    ## sum(c(-0.1,-0.1,+5)) > 0, so "2/3 agree" was graded not-reproduced.
    prof <- matrix(0, nrow = 3, ncol = 3,
                   dimnames = list(c("g1", "g2", "g3"), c("d1", "d2", "d3")))
    prof[1, ] <- c(-0.1, -0.1, 5)
    meta <- data.frame(dataset = c("d1", "d2", "d3"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(prof, meta, list(pw = c("g1", "g2", "g3")))
    expect_equal(out$majority_sign, -1L)
    expect_equal(out$concordance, 2 / 3)
    expect_equal(out$grade, "consistent")
})

test_that("datasets with no measured members do not inflate k", {
    prof <- matrix(0, nrow = 3, ncol = 3,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b", "c")))
    prof[1, ] <- c(-1, -1, -1)
    prof[, 3] <- NA_real_
    meta <- data.frame(dataset = c("a", "b", "c"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(prof, meta, list(pw = c("g1", "g2", "g3")))
    expect_equal(out$k, 2L)
    expect_equal(out$n_biosamples, 2L)
    expect_equal(out$grade, "consistent")   # k < 3, so reproduces cannot be assigned
})

test_that("finite member floors apply separately to every dataset", {
    prof <- matrix(c(-1, NA, Inf), 3, 3,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b", "c")))
    meta <- data.frame(dataset = c("a", "b", "c"), driver = "TF1",
                       biosample = c("A", "B", "C"))
    out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
    expect_equal(out$k, 0L)
    expect_equal(out$n_biosamples, 0L)
    expect_equal(out$grade, "not-computable")
    prof[, 1] <- -1
    out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
    expect_equal(out$k, 1L)
    expect_equal(out$grade, "not-computable")
})

test_that("missing biosamples never demonstrate independent reproduction", {
    prof <- matrix(-1, 3, 3,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b", "c")))
    meta <- data.frame(dataset = c("a", "b", "c"), driver = "TF1")
    out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
    expect_true(is.na(out$n_biosamples))
    expect_equal(out$grade, "consistent")
    for (bio in list(c("A", NA, NA), c("A", "", "  "), rep(NA_character_, 3))) {
        meta$biosample <- bio
        out <- mechanism_reproduction(prof, meta, list(pw = rownames(prof)))
        expect_lte(out$n_biosamples, 1L)
        expect_equal(out$grade, "consistent")
    }
})

test_that("dataset identities cannot inflate reproduction counts", {
    prof <- matrix(-1, 3, 2,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b")))
    meta <- data.frame(dataset = c("a", "b"), driver = "TF1", biosample = c("A", "B"))
    out <- mechanism_reproduction(prof, rbind(meta, meta[1, ]), list(pw = rownames(prof)))
    expect_equal(out$k, 2L)
    expect_equal(out$grade, "consistent")
    conflict <- meta[1, ]; conflict$biosample <- "C"
    expect_error(mechanism_reproduction(prof, rbind(meta, conflict), list(pw = rownames(prof))),
                 "inconsistent duplicate")
    colnames(prof) <- c("a", "a")
    expect_error(mechanism_reproduction(prof, meta, list(pw = rownames(prof))), "unique")
})

test_that("reproduction rows are matched through driver aliases", {
    ## GPSAdb returns NR3C1, but a user's table may use the alias GR. Without alias
    ## resolution the default wording silently comes back, and this is exactly the
    ## driver in two grade-B cases in P4 (airway / GSE286890).
    rep <- data.frame(driver = "GR", pathway = "pw", k = 3L,
                      n_biosamples = 3L, n_genes = 5L, concordance = 1,
                      majority_sign = -1L, grade = "reproduces",
                      stringsAsFactors = FALSE)
    hit <- .pe_reproduction_row("NR3C1", "pw", rep)
    expect_equal(nrow(hit), 1L)
    st <- .pe_single_line_statement("NR3C1", 5L, "D1", 0.01, hit)
    expect_match(st, "reproduces across 3 independent datasets", fixed = TRUE)
    ## The reverse direction must work too
    rep2 <- rep; rep2$driver <- "NR3C1"
    expect_equal(nrow(.pe_reproduction_row("GR", "pw", rep2)), 1L)
    expect_error(.pe_reproduction_row("NR3C1", "pw", data.frame(x = 1)),
                 "driver")
})
