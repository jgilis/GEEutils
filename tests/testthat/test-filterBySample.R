library(scuttle)

## Mock up data
set.seed(42)
n_cells <- 200L
n_genes <- 100L
ref <- mockSCE(ncells = n_cells, ngenes = n_genes)
ref$patient_id <- factor(rep(paste0("patient", 1:8), each = ncol(ref) / 8))

test_that("filterBySample() works", {
    n_drop <- 10L # number of genes that should be removed

    ## Set some genes to 0 for all but 1 patient
    counts(ref)[seq_len(n_drop), -which(ref$patient_id == "patient1")] <- 0
    to_drop <- rownames(ref)[seq_len(n_drop)]

    out <- filterBySample(ref, sample_id = "patient_id", size = 2)
    expect_equal(nrow(out), nrow(ref) - n_drop)
    expect_false(any(to_drop %in% rownames(out)))
})
