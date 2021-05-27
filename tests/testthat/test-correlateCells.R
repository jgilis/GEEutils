library(scuttle)

n_groups <- 3L
n_cells <- 50L

set.seed(42)
sce <- mockSCE(ncells = n_cells, ngenes = 100)
sce$subject <- sample(LETTERS[1:n_groups], n_cells, replace = TRUE)
sce <- logNormCounts(sce)

x <- logcounts(sce)

test_that("Vanilla correlateCells() works", {
    ## Hack around scran::correlatePairs to correlate cells instead of genes
    ref <- scran::correlatePairs(t(x))

    ## grouping = NULL, type = "within" calculates all pairwise correlations
    out <- correlateCells(sce, grouping = NULL, type = "within")
    expect_equal(sort(out$rho), sort(ref$rho))
    expect_named(out, c("cell1", "cell2", "rho", "p.value", "FDR"))
})

test_that("Within-group correlations work", {
    ## Calculate all pairwise correlations for each subject separately
    ref_A <- correlateCells(sce[, sce$subject == "A"])
    ref_B <- correlateCells(sce[, sce$subject == "B"])
    ref_C <- correlateCells(sce[, sce$subject == "C"])
    out <- correlateCells(sce, grouping = "subject", type = "within")
    expect_equal(out$rho, c(ref_A$rho, ref_B$rho, ref_C$rho))

    ## This only works for SCE inputs
    expect_error(correlateCells(x, grouping = "subject"))
})

test_that("Within-group correlations work for single-cell groups", {
    ## Use new subject for random single cell
    colData(sce)[sample(ncol(sce), 1), "subject"] <- "ZZZ"
    cells_per_subject <- table(sce$subject)
    n_pairs <- sum(vapply(cells_per_subject, choose, numeric(1), k = 2))
    out <- correlateCells(sce, grouping = "subject", type = "within")
    expect_equal(nrow(out), n_pairs)
})

test_that("Between-group correlations work", {
    ## Calculate how many total cell pairs we should get
    n_cells <- combn(table(sce$subject), 2)
    n_pairs <- sum(apply(n_cells, 2, prod))
    out <- correlateCells(sce, grouping = "subject", type = "between")
    expect_equal(nrow(out), n_pairs)
    expect_named(out, c("cell1", "cell2", "rho", "p.value", "FDR"))
})
