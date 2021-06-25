library(scuttle)

n_groups <- 3L
n_cells <- 50L

set.seed(42)
sce <- mockSCE(ncells = n_cells, ngenes = 100)
sce$subject <- sample(LETTERS[1:n_groups], n_cells, replace = TRUE)
sce <- logNormCounts(sce)

x <- logcounts(sce)

test_that("correlateCells() works with `grouping = NULL`", {
    ## Hack around scran::correlatePairs to correlate cells instead of genes
    ref <- scran::correlatePairs(t(x))

    ## grouping = NULL assumes all cells from single group and calculates all
    ## pairwise correlations
    out <- correlateCells(sce, grouping = NULL, n_pairs = NULL)
    expect_equal(sort(out$rho), sort(ref$rho))
    expect_named(out, c("cell1", "cell2", "rho", "p.value", "FDR", "type"))
})

test_that("correlateCells() works with both types", {
    ## Calculate how many total cell pairs we should get
    cells_per_subject <- table(sce$subject)
    within_pairs <- sum(vapply(cells_per_subject, choose, numeric(1), k = 2))
    cells_per_pair <- combn(cells_per_subject, 2)
    between_pairs <- sum(apply(cells_per_pair, 2, prod))

    out <- correlateCells(sce, grouping = "subject")
    expect_equal(nrow(out), within_pairs + between_pairs)
    expect_equal(sum(out$type == "within"), within_pairs)
    expect_equal(sum(out$type == "between"), between_pairs)
    expect_named(out, c("cell1", "cell2", "rho", "p.value", "FDR", "type"))
})

test_that("Within-group correlations work for single-cell groups", {
    ## Use new subject for random single cell
    colData(sce)[sample(ncol(sce), 1), "subject"] <- "ZZZ"
    cells_per_subject <- table(sce$subject)
    within_pairs <- sum(vapply(cells_per_subject, choose, numeric(1), k = 2))
    out <- correlateCells(sce, grouping = "subject")
    expect_equal(nrow(subset(out, type == "within")), within_pairs)
})

test_that("Subsetting cell pairs works", {
    n_pairs <- 10
    out <- correlateCells(sce, grouping = "subject", n_pairs = 10)
    expect_named(out, c("cell1", "cell2", "rho", "p.value", "FDR", "type"))
    ## Note: only the within-pairs will be exactly equal to n_pairs,
    ## for the between-pairs some deviation is possible due to the sampling
    ## strategy
    expect_equal(nrow(subset(out, type == "within")), n_pairs)
})
