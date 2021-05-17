library(scuttle)

n_groups <- 3L
n_cells <- 50L

set.seed(42)
sce <- mockSCE(ncells = n_cells, ngenes = 100)
sce$subject <- sample(LETTERS[1:n_groups], n_cells, replace = TRUE)
sce <- logNormCounts(sce)

x <- logcounts(sce)

## Hack around scran::correlatePairs to correlate cells instead of genes
ref <- scran::correlatePairs(t(x))

test_that("correlateCells() works", {
    ## grouping = NULL, type = "within" calculates all pairwise correlations
    out <- correlateCells(sce, grouping = NULL, type = "within")
    expect_equal(sort(out), sort(ref$rho))
})

test_that("Within-group correlations work", {
    ref_A <- correlateCells(sce[, sce$subject == "A"])
    ref_B <- correlateCells(sce[, sce$subject == "B"])
    ref_C <- correlateCells(sce[, sce$subject == "C"])
    out <- correlateCells(sce, grouping = "subject", type = "within")
    expect_equal(out, c(ref_A, ref_B, ref_C))

    expect_error(correlateCells(x, grouping = "subject"))
})
