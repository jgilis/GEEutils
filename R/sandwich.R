# TODO: add some more details regarding "clustered covariance matrices" and within-subject correlation
# FIXME: shouldn't this function always fail when number of parameters is lower than number of subjects?

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
#'   estimator is used (equivalent to [sandwich::vcovHC()]).
#' @param type Character string specifying which type of small-sample bias
#'   adjustment to be applied. Possible values are `"HC0"` to `"HC3"` or
#'   `"LiRedden"`. For details on the `"HC*"` adjustments see
#'   [sandwich::vcovCL()], for `"LiRedden"` see "Details" below. The default is
#'   to use `"LiRedden"`.
#' @param cadjust Logical, whether to apply cluster bias adjustment. Where
#'   "cluster" refers to the subject-level grouping of the observations. Only
#'   relevant if `subject_id` is provided. See [sandwich::vcovCL()] for details.
#'   Default: `TRUE`. When using `type = "LiRedden"`, this will be set to
#'   `FALSE` internally.
#' @param fix Logical, whether to fix the covariance matrix to be positive
#'   semi-definite in case it's not. See [sandwich::vcovCL()] for details.
#'   Default: `FALSE`. Only relevant for `type = "HC2"` or `"HC3"`.
#' @param use_T Logical, should a \eqn{t}-test be used to calculate the
#'   p-values? If `FALSE` (the default), will use a \eqn{z}-test.
#'
#' @details
#' ## Small sample size adjustments
#'
#' Note that 'sample size' here refers to the number of subjects, *not* the
#' number of cells.
#'
#' The `"LiRedden"` adjustment uses the degrees-of-freedom correction proposed
#' in Li & Redden (2015). Briefly, the variances of \eqn{\hat{\beta}} are
#' multiplied with \eqn{K/(K-p)}, where \eqn{K} is the number of subjects and
#' \eqn{p} is the number of regression parameters.
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
#' Li P. and Redden D. T. 'Small sample performance of bias-corrected sandwich
#' estimators for cluster-randomized trials with binary outcomes', *Statistics
#' in Medicine* __34__(2), 281–296 (2015).
#' doi: [10.1002/sim.6344](https://doi.org/10.1002/sim.6344).
#'
#' Zeileis A., Köll S. and Graham N. (2020). 'Various Versatile Variances: An
#' Object-Oriented Implementation of Clustered Covariances in R', *Journal of
#' Statistical Software*, __95__(1), 1–36.
#' doi: [10.18637/jss.v095.i01](https://doi.org/10.18637/jss.v095.i01).
#'
#' @examples
#' ## Mock up data set
#' library(scuttle)
#' sce <- mockSCE(ncells = 1000, ngenes = 500)
#'
#' ## Simulate 8 samples in two treatment groups, using an unpaired design
#' sce$Subjects <- gl(8, 125)
#' sce$Treatment <- gl(2, 4)[sce$Subjects]
#'
#' colData(sce)
#'
#' ## Fit model using available colData columns
#' model_fits <- fitGLM(sce, formula = ~Treatment)
#'
#' ## Test for Treatment effect using a clustered sandwich estimator with
#' ## Li-Redden adjustment
#' res <- glmSandwichTest(
#'     model_fits,
#'     coef = 2,
#'     subject_id = "Subjects",
#'     type = "LiRedden",
#'     use_T = TRUE  # use Wald t-tests
#' )
#' head(res$table)
#'
#' @export
#' @importFrom stats df.residual p.adjust pnorm pt
glmSandwichTest <- function(models, subject_id,
                            coef = NULL, contrast = NULL,
                            type = c("LiRedden", "HC0", "HC1", "HC2", "HC3"),
                            cadjust = FALSE, fix = FALSE,
                            use_T = FALSE) {

    type <- match.arg(type)

    ## If single model supplied, create list internally
    if (!is(models, "list")) {
        if (!is(models, "glm")) {
            stop(
                "`models` should be a single 'glm' or a list of 'glm' objects.",
                call. = FALSE
            )
        }
        models <- list(models)
    }

    d <- models[[1]]$data
    if (!is.null(subject_id) && !(subject_id %in% colnames(d))) {
        stop(sprintf("subject_id = '%s' not found.", subject_id), call. = FALSE)
    }

    coefficients <- .get_coefs(models = models)
    comparison <- .handle_comparison(
        coefficients = coefficients, coef = coef,
        contrast = contrast
    )
    effect_size <- comparison$effect_size
    comp_name <- comparison$comp_name

    sandwich_out <- .get_beta_vars(
        models = models, coef = coef, contrast = contrast,
        subject_id = subject_id, type = type, cadjust = cadjust, fix = fix
    )
    sandwich_var <- sandwich_out$beta_vars

    se <- sqrt(sandwich_var)
    wald_stats <- effect_size / se

    if (use_T) {
        m <- models[[1]]  # same for all fits
        if (!is.null(subject_id)) {
            df <- .get_between_subject_df(m, subject_id = subject_id)
        } else {
            df <- m$df.residual
        }
        pvals <- 2 * pt(abs(wald_stats), df = df, lower.tail = FALSE)
    } else {
        pvals <- 2 * pnorm(abs(wald_stats), lower.tail = FALSE)
    }
    fdr <- p.adjust(pvals, method = "fdr")

    gene_names <- names(models)
    tab <- data.frame(
        logFC = effect_size / log(2),
        SE = se,
        Wald = wald_stats,
        PValue = pvals,
        FDR = fdr,
        comparison = comp_name,
        row.names = gene_names
    )
    params <- c(
        subject_id = subject_id,
        sandwich_out$params,
        use_T = use_T
    )
    if (use_T) params$df <- df

    ## Return params and table
    list(params = params, table = tab)
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
    adjust_LR <- FALSE
    if (type == "LiRedden") {
        type <- "HC0"
        if (is.null(subject_id)) {
            warning(
                '`type = "LiRedden"` is only relevant when `subject_id` is provided.',
                '\nFalling back to `type = "HC0"` without further adjustment.'
            )
        } else {
            cadjust <- fix <- FALSE
            adjust_LR <- TRUE
            lr_adjustment <- .LR_adjustment(models[[1]], subject_id = subject_id)
        }
    }

    if (is.character(subject_id)) {
        stopifnot(length(subject_id) == 1)
        subject_id <- as.formula(paste("~", subject_id))
    }

    beta_vars <- vapply(models, FUN = .glm_sandwich_var,
        contrast = contrast, coef = coef,
        subject_id = subject_id, type = type,
        cadjust = cadjust, fix = fix,
        FUN.VALUE = numeric(1)  # expects just one coefficient or contrast!
    )

    if (adjust_LR) {
        type <- lr_adjustment$type
        beta_vars <- lr_adjustment$out * beta_vars
    }

    list(
        beta_vars = beta_vars,
        params = list(type = type, cadjust = cadjust, fix = fix)
    )
}


## Helper to compute sandwich variances of coefficient estimates for single model
## Returns diagonal elements of the covariance matrix V(beta)
#' @importFrom sandwich vcovCL
.glm_sandwich_var <- function(m, coef, contrast, subject_id, type, cadjust, fix) {
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


.LR_adjustment <- function(m, subject_id) {
    d <- m$data
    if (!(subject_id %in% names(d))) {
        stop(
            "`subject_id` '", subject_id, "' not found in fit data.",
            call. = FALSE
        )
    }
    K <- nlevels(factor(d[[subject_id]]))

    ## Do not count within-subject params
    n_within <- .count_within_params(m, subject_id = subject_id)
    p <- length(m$coefficients) - n_within

    if (K <= p) {
        warning(
            "Number of subjects (", K, "), ",
            "lower than number of parameters (", p, ").",
            "\nLi-Redden adjustment not possible.",
            "\nFalling back to `type = 'HC0'`.",
            call. = FALSE
        )
        type <- "HC0"
        out <- 1 # no LR adjustment == HC0
    } else {
        type <- "LiRedden"
        out <- K / (K - p)
    }
    list(out = out, type = type)
}


## Helpers to compute the between-subject degrees of freedom
.get_between_subject_df <- function(m, subject_id) {
    d <- m$data

    ## df = nr. of subjects - nr. of between-subject parameters
    K <- length(unique(d[[subject_id]]))
    p <- length(m$coefficients)
    n_within_params <- .count_within_params(m, subject_id = subject_id)
    df <- K - (p - n_within_params)

    if (df <= 0) {
        stop("Non-positive degrees of freedom.\n",
            "Too many between-subject parameters (", p - n_within_params, ") ",
            "for number of subjects (", K, ").",
            call. = FALSE
        )
    }

    df
}

#' @importFrom stats formula
.count_within_params <- function(m, subject_id) {
    ## Get variables from model formula
    ff <- formula(m)
    all_vars <- all.vars(ff)

    ## Remove 'y' and 'offsets' variables
    ## Don't need to consider those
    all_vars <- all_vars[!(all_vars %in% c("y", "offsets"))]

    ## Check which are within-variables
    within_vars <- .get_within_vars(
        all_vars,
        .data = m$data, subject_id = subject_id
    )

    if (!length(within_vars)) {
        return(0)
    }

    ## Count parameters in which the within_vars appear, also the interactions
    ## These would all be considered within-subject variables
    coefs <- names(m$coefficients)

    occurrences <- vapply(within_vars, grepl, logical(length(coefs)), x = coefs)
    occurrences <- apply(occurrences, 1, any)  # avoid double counting
    sum(occurrences)
}

.is_within_var <- function(var, .data, subject_id) {
    any(rowSums(table(.data[[subject_id]], .data[[var]]) > 0) > 1)
}

.get_within_vars <- function(vars, .data, subject_id) {
    is_within <- vapply(vars, .is_within_var, logical(1),
        .data = .data, subject_id = subject_id
    )
    vars[is_within]
}
