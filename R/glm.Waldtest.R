# Extract model coefficients from a glm
.glm_getCoef <- function(model) {
    model$coefficients
}

# Compute estimates for target contrast
.glm_getEstimates <- function(model, contrast) {
    coef <- .glm_getCoef(model)
    if (is.null(coef)) {
        coef <- rep(NA, times = nrow(contrast))
    }
    contrast %*% coef
}

# Compute the standard error on a target model estimate
#' @importFrom methods is
.glm_varContrast <- function(model, contrast) {
    if (!any(is(model, "fitError"))) {
        sx <- summary(model) # maybe I can also only compute vcov to gain speed
        if (nrow(sx$cov.scaled) == length(contrast)) {
            return(sqrt(diag(t(contrast) %*% sx$cov.scaled %*% contrast))) # scaled by dispersion
        }
    }
    NA
}

# Compute the robust variance on a target model estimate
#' @importFrom methods is
#' @importFrom sandwich sandwich
.glm_sandwichVarContrast <- function(model, contrast, adjust) {
    if (!any(is(model, "fitError"))) {
        # if adjust == TRUE: small sample adjustment, multiply by n/(n-k)
        sanwich_var <- sandwich(model, adjust = adjust)
        if (nrow(sanwich_var) == length(contrast)) {
            return(diag(t(contrast) %*% sanwich_var %*% contrast))
        }
    }
    NA
}

# Extract the degrees of freedom from a glm
.glm_getDf <- function(model) {
    model$df.residual
}


#' Perform a Wald test on a glm with the possible of adopting robust standard
#' errors
#'
#' @description Perform a Wald test on a glm with the possible of adopting
#'   robust standard errors
#'
#' @param models A list of generalized linear model. If a single model, make
#'   sure to input it as a list by `list(glm)`.
#'
#' @param contrast A `matrix` specifying the contrast, i.e. combinations of
#'   model parameters, of interest.
#'
#' @param sandwich A `logical` variable; should robust standard errors be
#'   computed (sandwich procedure by Liang and Zeger). Default is TRUE.
#'
#' @param adjust A `logical` variable; should a finite sample adjustment be
#'   made? This amounts to multiplication with n/(n-k) where n is the number of
#'   observations and k the number of estimated parameters.
#'
#' @return A `Dataframe` containing the requested model parameter estimates,
#'   (robust) standard errors, degrees of freedom, Wald test statistics and
#'   p-values.
#'
#' @rdname glm.Waldtest
#'
#' @author Jeroen Gilis
#'
#' @examples
#' ## Simulate single gene across 10 patients in 2 groups
#' n_cells <- 500
#' gene <- rpois(n_cells, 3)
#' group_id <- factor(rep(c("control", "treat"), each = n_cells / 2))
#' patient_id <- factor(rep(paste0("patient", 1:10), each = n_cells / 10))
#' data <- data.frame(gene, group_id, patient_id)
#'
#' ## Set up design and contrast
#' design <- model.matrix(~group_id, data = data)
#' L <- matrix(0, ncol = 1, nrow = ncol(design))
#' rownames(L) <- colnames(design)
#' L["group_idtreat", 1] <- 1
#'
#' ## Fit Poisson model
#' pois_model <- glm(gene ~ group_id, family = poisson, data = data)
#'
#' ## Wald test with sandwich estimator and small sample-size adjustment
#' res <- glm.Waldtest(
#'     models = list(pois_model),
#'     contrast = L[, 1],
#'     sandwich = TRUE,
#'     adjust = TRUE
#' )
#'
#' @export
glm.Waldtest <- function(models, contrast, sandwich = TRUE, adjust = FALSE) {
    # FIXME: replace `sapply` with `vapply`
    # cfr. https://bioconductor.org/developers/package-guidelines/#rcode
    # TODO: convert 'models' to list internally? Makes using a single model
    # bit more intuitive
    estimates <- sapply(models, .glm_getEstimates, contrast = contrast)
    if (sandwich) {
        var <- sapply(models, .glm_sandwichVarContrast, contrast = contrast, adjust) # uses sandwich SEs
    } else {
        var <- sapply(models, .glm_varContrast, contrast = contrast) # uses GLM SEs
    }
    dfs <- sapply(models, .glm_getDf)
    
    se <- sqrt(var) # standard error
    W_stats <- estimates / se # Wald statistics
    pvals <- 2 * pt(abs(W_stats), df = dfs, lower.tail = FALSE)

    data.frame(estimates, se, dfs, W_stats, pvals)
}
