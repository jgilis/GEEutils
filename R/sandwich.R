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

    # TODO: handle type == "LiRedden"

    sandwich_var <- vapply(models, FUN = .glm_sandwich_var,
        contrast = contrast, coef = coef,
        subject_id = subject_id, type = type,
        cadjust = cadjust, fix = fix,
        FUN.VALUE = numeric(NCOL(effect_size))
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
