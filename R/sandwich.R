# TODO: add some more details regarding "clustered covariance matrices" and within-subject correlation
# TODO: add details regarding Li-Redden adjustment
# TODO: replace example with more sensible case

#' GLM Wald tests using sandwich estimators
#'
#' Perform Wald tests GLM coefficients using sandwich covariance estimators.
#' When `subject_id` is specified, will use the __clustered covariance
#' estimator__, implemented in [sandwich::vcovCL()], to account for
#' within-subject correlations.
#'
#' @param models List of `glm` models, resulting from [fitGLM()].
#' @param coef Character string or numeric scalar specifying the coefficient to
#'   be tested. Ignored if `contrast` is provided.
#' @param contrast Integer vector specifying a linear combination of
#'   coefficients to be tested. Should have length equal to the number of
#'   coefficients in the model. Takes precedence over `coef`.
#' @param subject_id Optional character string specifying a column in the
#'   original data indicating subject IDs. When this is specified, observations
#'   are considered to be correlated within each subject and
#'   [sandwich::vcovCL()] is used to estimate the covariance matrix. When
#'   `subject_id = NULL` (the default), a heteroskedastic-consistent sandwich
#'   estimator is used (equivalent to [sandwich::vcovHC()]). Can also provide a
#'   [formula][stats::formula].
#' @param type Character string specifying which type of small-sample bias
#'   adjustment to be applied. Possible values are `"HC0"` to `"HC3"` or
#'   `"LiRedden"`. For details on the `"HC*"` adjustments see
#'   [sandwich::vcovCL()], for `"LiRedden"` see "Details" below. The default is
#'   to use `"HC0"`, which corresponds to White's estimator.
#' @param cadjust Logical, whether to apply cluster bias adjustment. Where
#'   "cluster" refers to the subject-level grouping of the observations. Only
#'   relevant if `subject_id` is provided. See [sandwich::vcovCL()] for details.
#'   Default: `TRUE`.
#' @param fix Logical, whether to fix the covariance matrix to be positive
#'   semi-definite in case it's not. See [sandwich::vcovCL()] for details.
#'   Default: `FALSE`.
#'
#' @details
#' The `"LiRedden"` adjustment......................
#'
#' @return
#' A `data.frame` with a row for each gene and the following columns:
#'
#' * `logFC`: log2 Fold-change estimates for the given comparison
#' * `SE`: standard errors from the requested sandwich estimator type
#' * `Wald`: Wald statistics
#' * `PValue`
#' * `FDR`: FDR-adjusted p-values using Benjamini-Hochberg's method
#' * `comparison`: the coefficient or contrast tested
#'
#' @seealso [sandwich::vcovCL()]
#' @references
#' Zeileis (2020)
#' Li & Redden (2015)
#'
#' @examples
#' ## Mock up data set
#' library(scuttle)
#' sce <- mockSCE(ncells = 100, ngenes = 500)
#' colData(sce)
#'
#' ## Fit model using available colData columns
#' model_fits <- fitGLM(sce, ~ Mutation_Status + Treatment)
#'
#' ## Test for Treatment effect using a HC sandwich estimator with HC3 adjustment
#' res_HC3 <- glmSandwichTest(
#'     model_fits,
#'     coef = "Treatmenttreat2",
#'     type = "HC3"
#' )
#' head(res_HC3)
#'
#' ## Test for Treatment effect using a Clustered sandwich estimator with HC3 adjustment
#' ## Assuming cells are correlated within each "Cell_Cycle" level
#' res_CL <- glmSandwichTest(
#'     model_fits,
#'     coef = "Treatmenttreat2",
#'     subject_id = "Cell_Cycle", type = "HC3"
#' )
#' head(res_CL)
#'
#' @export
#' @importFrom stats df.residual p.adjust pnorm
glmSandwichTest <- function(models, coef = NULL, contrast = NULL,
                            subject_id = NULL,
                            type = c("HC0", "HC1", "HC2", "HC3", "LiRedden"),
                            cadjust = TRUE, fix = FALSE) {

    type <- match.arg(type)

    ## If single model supplied, create list internally
    if (is(models, "glm")) {
        models <- list(models)
    }

    coefficients <- .get_coefs(models = models)
    comparison <- .handle_comparison(
        coefficients = coefficients, coef = coef,
        contrast = contrast
    )
    effect_size <- comparison$effect_size
    comp_name <- comparison$comp_name

    sandwich_var <- .get_beta_vars(
        models = models, coef = coef, contrast = contrast,
        subject_id = subject_id, type = type, cadjust = cadjust, fix = fix
    )

    se <- sqrt(sandwich_var)
    wald_stats <- effect_size / se

    pvals <- 2 * pnorm(abs(wald_stats), lower.tail = FALSE)
    fdr <- p.adjust(pvals, method = "fdr")

    gene_names <- names(models)
    data.frame(
        logFC = effect_size / log(2),
        SE = se,
        Wald = wald_stats,
        PValue = pvals,
        FDR = fdr,
        comparison = comp_name,
        row.names = gene_names
    )
}


#' @importFrom stats coef
.get_coefs <- function(models) {
    n_coef <- length(coef(models[[1]]))
    t(vapply(models, coef, numeric(n_coef)))
}


## Helper to compute effect size of `coef` or `contrast`
.handle_comparison <- function(coefficients, coef, contrast) {
    coef_names <- colnames(coefficients)

    if (is.null(contrast)) {
        if (is.null(coef)) {
            stop("Need to provide either `coef` or `contrast`.", call. = FALSE)
        } else {
            stopifnot(length(coef) == 1) # currently only support testing single coef

            if (is.character(coef)) {
                coef_present <- coef %in% coef_names
                if (!all(coef_present)) {
                    stop("The following `coef` names are missing:",
                         "\n\t", paste(coef[!coef_present], collapse = ", "),
                         call. = FALSE
                    )
                }
                comparison <- coef
                coef <- match(coef, coef_names)
            } else if (is.numeric(coef)) {
                stopifnot(coef <= length(coefficients))
                comparison <- coef_names[coef]
            } else {
                stop("`coef` must be either numeric or character, not ",
                     typeof(coef),
                     call. = FALSE
                )
            }
            out <- coefficients[, coef]
        }
    } else {
        contrast <- as.matrix(contrast)

        stopifnot(ncol(contrast) == 1) # currently only support testing single contrast

        if (nrow(contrast) != ncol(coefficients)) {
            stop("`nrow(contrast)` should be equal to the number of coefficients.",
                 call. = FALSE
            )
        }

        ncontrasts <- qr(contrast)$rank
        if (ncontrasts == 0) stop("All contrasts are zero.", call. = FALSE)
        out <- drop(coefficients %*% contrast)

        if (is.null(colnames(contrast))) {
            comparison <- vapply(seq_len(ncol(contrast)), function(i) {
                contr_vec <- contrast[, i]
                idx <- which(contr_vec != 0)
                comparison <- paste(contr_vec[idx], coef_names[idx], sep = "*")
                comparison <- gsub("^([-]?)1\\*", "\\1", comparison)
                paste(comparison, collapse = " + ")
            }, character(1))
        } else {
            comparison <- colnames(contrast)
        }
    }
    list(effect_size = out, comp_name = comparison)
}


.get_beta_vars <- function(models, coef, contrast, subject_id, type, cadjust, fix) {
    # TODO: handle type == "LiRedden"
    out <- vapply(models, FUN = .glm_sandwich_var,
        contrast = contrast, coef = coef,
        subject_id = subject_id, type = type,
        cadjust = cadjust, fix = fix,
        FUN.VALUE = numeric(1)  # expects just one coefficient or contrast!
    )

    out
}


## Helper to compute sandwich variances of coefficient estimates for single model
## Returns diagonal elements of the covariance matrix V(beta)
#' @importFrom sandwich vcovCL
.glm_sandwich_var <- function(m, coef, contrast, subject_id, type, cadjust, fix) {
    if (is.character(subject_id)) {
        stopifnot(length(subject_id) == 1)
        subject_id <- as.formula(paste("~", subject_id))
    }
    v <- vcovCL(m,
        cluster = subject_id, sandwich = TRUE, type = type,
        cadjust = cadjust, fix = fix
    )
    if (!is.null(contrast)) {
        v <- crossprod(contrast, v %*% contrast)
        return(diag(v))
    } else if (!is.null(coef)) {
        ## Only return for requested coefficients
        return(diag(v)[coef])
    }
    diag(v)
}
