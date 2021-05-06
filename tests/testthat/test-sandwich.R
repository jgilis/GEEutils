## Mock up data
set.seed(42)
n_cells <- 50L
n_genes <- 100L
n_subjects <- 5L

sce <- scuttle::mockSCE(ncells = n_cells, ngenes = n_genes)
sce$subject_id <- rep(paste0("subject", 1:n_subjects), each = ncol(sce) / n_subjects)
model_fits <- fitGLM(sce, ~ Mutation_Status + Treatment)

test_that("Vanilla glmSandwichTest() works", {
    ## Uses 'Li-Redden' adjustment by default
    out <- glmSandwichTest(model_fits, subject_id = "subject_id", coef = 3)

    pars <- out$params
    expect_identical(pars$type, "LiRedden")
    expect_identical(pars$subject_id, "subject_id")
    expect_false(pars$cadjust)
    expect_false(pars$fix)

    tab <- out$table
    expect_s3_class(tab, "data.frame")
    expect_identical(nrow(tab), length(model_fits))
})

test_that("glmSandwichTest() works with single model fit", {
    m <- model_fits[[1]]
    out <- glmSandwichTest(m, subject_id = "subject_id", coef = 3)
    tab <- out$table
    expect_s3_class(tab, "data.frame")
    expect_identical(nrow(tab), 1L)
})

test_that("glmSandwichTest() works with contrasts", {
    contr <- c(-1, 0, 1)
    out <- glmSandwichTest(model_fits, subject_id = "subject_id", contrast = contr)
    tab <- out$table
    expect_s3_class(tab, "data.frame")
    expect_identical(nrow(tab), length(model_fits))
})

test_that("glmSandwichTest() fails and warns when necessary", {
    expect_error(
        glmSandwichTest(model_fits, subject_id = "subject_id", coef = "random")
    )
    expect_error(
        glmSandwichTest(model_fits, subject_id = "subject_id", coef = 1e5)
    )
    expect_error(
        glmSandwichTest(model_fits, subject_id = "subject_id", coef = TRUE)
    )
    contr <-  cbind(1:3, 4:6)
    expect_error(
        glmSandwichTest(model_fits, subject_id = "subject_id", contrast = contr)
    )
    expect_error(glmSandwichTest(model_fits, subject_id = "subject_id"))
    expect_error(glmSandwichTest(model_fits, subject_id = "bla", coef = 3))
    expect_warning(
        out <- glmSandwichTest(
            model_fits, subject_id = NULL, coef = 3, type = "LiRedden"
        )
    )
    expect_identical(out$params$type, "HC0")
})
