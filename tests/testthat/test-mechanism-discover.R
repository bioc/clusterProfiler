## Decision rules and guards for the discovery gate (mechanism_discover).
##
## This layer touches no network and no external table, so all of it is testable.

.make_net <- function(df) data.frame(source = df[[1]], target = df[[2]],
                                     weight = df[[3]], stringsAsFactors = FALSE)

test_that("a factor whose targets moved ranks first", {
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(c(rep(1.2, 6), rep(0, 34)), genes)
    ## Spell out data.frame: `.make_net(list(...))` recycles on unequal vector
    ## lengths and interleaves the two TFs' targets into one set (the first version
    ## of this test failed exactly that way).
    net <- data.frame(source = rep(c("TF1", "TF2"), each = 6),
                      target = c(genes[1:6], genes[7:12]), weight = 1,
                      stringsAsFactors = FALSE)
    out <- mechanism_discover(signed, net, top = 2)
    expect_equal(out$source[1], "TF1")
    expect_equal(out$direction[1], "activates")
    expect_gt(out$score[1], out$score[2])
})

test_that("self-targets are dropped, and that is what keeps knockouts honest", {
    ## In a knockout experiment the TF's own log2FC is the strongest; a self-loop
    ## pushes its score to the top. This is not theoretical: 62% of the TFs in
    ## CollecTRI carry a self-loop, and removing them is worth 5 points.
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(c(-3, rep(1.2, 6), rep(0, 33)), genes)
    ## In real data a self-loop is "the TF lists itself as a target" - so the TF name
    ## must be one of the gene names. g1 is TF1: its targets = itself (g1) + g2..g7;
    ## TF2's targets = g2..g7.
    net <- data.frame(
        source = c("g1", rep("g1", 6), rep("g2", 6)),
        target = c("g1", genes[2:7], genes[2:7]),
        weight = 1, stringsAsFactors = FALSE)
    expect_message(out <- mechanism_discover(signed, net, top = 2), "self-target")
    ## With self-loops removed TF1 and TF2 have identical target sets and tied
    ## scores; the point is that TF1 no longer leads on its own log2FC
    expect_equal(out$n_targets[out$source == "g1"], 6L)
    keep <- suppressMessages(mechanism_discover(signed, net, top = 2, keep_self_targets = TRUE))
    expect_equal(keep$n_targets[keep$source == "g1"], 7L)
})

test_that("weights carry the direction of regulation", {
    genes <- paste0("g", 1:40)
    ## The targets go down overall, but the weights say "represses target" (-1)
    ## -> the factor should be called activating
    signed <- stats::setNames(c(rep(-1.2, 6), rep(0, 34)), genes)
    net <- .make_net(list("TF1", genes[1:6], rep(-1, 6)))
    out <- mechanism_discover(signed, net, top = 1)
    expect_equal(out$direction[1], "activates")
    expect_gt(out$score[1], 0)
})

test_that("factors with too few measured targets are skipped", {
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(stats::rnorm(40), genes)
    net <- data.frame(source = c(rep("TF1", 3), rep("TF2", 17)),
                      target = c(genes[1:3], genes[4:20]), weight = 1,
                      stringsAsFactors = FALSE)
    out <- mechanism_discover(signed, net, top = 5, min_targets = 5)
    expect_equal(nrow(out), 1)
    expect_equal(out$source, "TF2")
})

test_that("targets outside the measured genes are dropped, not zero-filled", {
    ## Zero-filling treats "not measured" as "unchanged" and dilutes the score.
    ## It must be computed on the measured intersection.
    genes <- paste0("g", 1:40)
    ## Only 5 genes are covered by the network, so they must differ among themselves,
    ## otherwise the variance is 0 (this guard is correct)
    signed <- stats::setNames(c(1, 1.1, 0.9, 1.05, 0.95, rep(0, 35)), genes)
    net <- data.frame(source = "TF1",
                      target = c(genes[1:5], "not_measured_1", "not_measured_2"),
                      weight = 1, stringsAsFactors = FALSE)
    out <- mechanism_discover(signed, net, top = 1)
    expect_equal(out$n_targets, 5)
    expect_equal(attr(out, "n_genes"), 40)
})

test_that("a network with no overlap fails loudly", {
    signed <- stats::setNames(c(a = 1, b = 2), c("a", "b"))
    net <- .make_net(list("TF1", c("x", "y", "z"), rep(1, 3)))
    expect_error(mechanism_discover(signed, net), "None of the network targets")
})

test_that("the shortlist is capped and ordered", {
    genes <- paste0("g", 1:60)
    signed <- stats::setNames(stats::rnorm(60), genes)
    net <- .make_net(list(rep(paste0("TF", 1:5), each = 8),
                          rep(genes[1:40], 5), rep(1, 40)))
    out <- mechanism_discover(signed, net, top = 3)
    expect_equal(nrow(out), 3)
    expect_equal(out$rank, 1:3)
    expect_true(all(diff(out$score) <= 0))
})

test_that("rank_by = abs puts the strongest negative score first", {
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(c(rep(-1.2, 6), rep(0, 34)), genes)
    net <- data.frame(source = c(rep("TF_down", 6), rep("TF_up", 6)),
                      target = c(genes[1:6], genes[7:12]),
                      weight = c(rep(1, 6), rep(-1, 6)),
                      stringsAsFactors = FALSE)
    by_score <- mechanism_discover(signed, net, top = 2, rank_by = "score")
    by_abs <- mechanism_discover(signed, net, top = 2, rank_by = "abs")
    expect_equal(by_score$source[1], "TF_up")
    expect_equal(by_abs$source[1], "TF_down")
})

test_that("share uses the whole measured vector and tested counts are honest", {
    genes <- paste0("g", 1:40)
    signed <- stats::setNames(stats::rnorm(40), genes)
    net <- data.frame(source = rep(c("TF1", "TF2", "TF3"), c(6, 6, 2)),
                      target = c(genes[1:6], genes[7:12], genes[13:14]),
                      weight = 1, stringsAsFactors = FALSE)
    out <- mechanism_discover(signed, net, top = 5, min_targets = 5)
    ## TF3 has only 2 target genes and is skipped by min_targets: tested should be 2, not 3
    expect_equal(attr(out, "n_factors_tested"), 2L)
    expect_equal(attr(out, "n_factors_in_network"), 3L)
    expect_equal(attr(out, "n_genes"), 40L)
    expect_equal(out$share[out$source == "TF1"], 6 / 40)
})

test_that("duplicated gene names are collapsed by largest |log2FC|", {
    genes <- c(paste0("g", 1:39), "g1")
    signed <- stats::setNames(c(2, rep(0, 39)), genes)
    signed[40] <- 0.1
    net <- data.frame(source = "TF1", target = paste0("g", 1:6), weight = 1,
                      stringsAsFactors = FALSE)
    expect_message(out <- mechanism_discover(signed, net, top = 1), "duplicated")
    expect_equal(attr(out, "n_genes"), 39L)
})

test_that("conflicting duplicate weights cannot be resolved by row order", {
    signed <- setNames(c(rep(1, 5), rep(0, 15)), paste0("g", 1:20))
    net <- data.frame(source = "TF", target = c(paste0("g", 1:5), "g1"),
                      weight = c(rep(1, 5), -1))
    expect_error(mechanism_discover(signed, net), "Conflicting weights")
    expect_error(mechanism_discover(signed, net[nrow(net):1, ]), "Conflicting weights")
    clean <- net[1:5, ]
    expect_equal(mechanism_discover(signed, rbind(clean, clean)),
                 mechanism_discover(signed, clean))
})

test_that("non-finite measurements do not contaminate the correlation", {
    signed <- setNames(c(1, 2, 3, 4, 5, rep(0, 14), Inf), paste0("g", 1:20))
    net <- data.frame(source = "TF", target = paste0("g", 1:5), weight = 1)
    expect_message(out <- mechanism_discover(signed, net), "non-finite")
    expect_equal(attr(out, "n_genes"), 19L)
    expect_true(is.finite(out$score))
})
