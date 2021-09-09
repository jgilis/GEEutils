## Mock up data
set.seed(42)
n_cells <- 80L
n_genes <- 100L
n_subjects <- 5L

sce <- scuttle::mockSCE(ncells = n_cells, ngenes = n_genes)
sce$subject_id <- gl(8, 10)
sce$Treatment <- gl(2, 4)[sce$subject_id]
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


## Replace subject_id column with only 2 subjects, which is lower than the 3
## parameters in the GLM fits
n_subjects <- 2
sce2 <- sce
sce2$subject_id <- rep(paste0("subject", 1:n_subjects), each = ncol(sce) / n_subjects)
model_fits2 <- fitGLM(sce2, ~ Mutation_Status + Treatment)

## Reference: HC0 adjustment, if K < p, LiRedden should be identical to this
ref_HC0 <- glmSandwichTest(
    model_fits2, subject_id = "subject_id", coef = 3, type = "HC0"
)

# FIXME: shouldn't this function always fail when number of parameters is lower than number of subjects?
test_that("LiRedden adjustment handles K <= p correctly", {
    expect_warning(
        out <- glmSandwichTest(
            model_fits2,
            subject_id = "subject_id", coef = 3, type = "LiRedden"
        )
    )
    expect_identical(out$params$type, "HC0")
    expect_identical(out, ref_HC0)
})


# t-tests -----------------------------------------------------------------

## Refit models but now with 2 within-subject variables and an interaction
model_fits3 <- fitGLM(sce, ~ Mutation_Status + Cell_Cycle * Treatment)

test_that("glmSandwichTest() works with t-tests", {
    out <- glmSandwichTest(model_fits3,
        subject_id = "subject_id",
        type = "LiRedden",
        coef = "Treatment2", use_T = TRUE
    )

    pars <- out$params
    expect_identical(pars$type, "LiRedden")
    expect_identical(pars$subject_id, "subject_id")

    expect_true(pars$use_T)

    ## df = 6 = 8 subjects - 2 between-subject variables (intercept + Treatment)
    expect_identical(pars$df, 6L)
})

## Case with too many parameters
sce3 <- sce
trtB <- rep(c("yes", "no"), length.out = nlevels(sce3$subject_id))
sce3$TreatmentB <- trtB[sce$subject_id]
trtC <- rep(c("yes", "no"), length.out = nlevels(sce3$subject_id))
sce3$TreatmentC <- trtC[sce$subject_id]

model_fits4 <- fitGLM(sce3, ~ Treatment * TreatmentB * TreatmentC)

test_that("Sandwich t-test fails with too few subjects", {
    expect_error(
        glmSandwichTest(model_fits4,
            subject_id = "subject_id",
            type = "HC0",  # Use HC0 to avoid Li-Redden warning
            coef = "Treatment2", use_T = TRUE
        ),
        "Non-positive degrees of freedom."
    )
})


## Now with an interaction between two within-subject variabels
model_fits5 <- fitGLM(sce, ~ Mutation_Status * Cell_Cycle + Treatment)

test_that("Sandwich t-test handles within-subject interaction effects", {
    out <- glmSandwichTest(model_fits5,
        subject_id = "subject_id",
        type = "LiRedden",
        coef = "Treatment2", use_T = TRUE
    )

    pars <- out$params
    expect_true(pars$use_T)

    ## df should still be 6
    expect_identical(pars$df, 6L)
})
