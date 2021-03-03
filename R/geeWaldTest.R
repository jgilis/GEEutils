.gee_getCoef <- function(model) {
    if (class(model)[1] == "try-error") {
        return(NULL)
    }
    return(model$coefficients)
}

.gee_getEstimates <- function(model, contrast) {
    coef <- .gee_getCoef(model)
    if (is.null(coef)) {
        coef <- rep(NA, times = length(contrast))
    }
    contrast %*% coef
}

.gee_varContrast <- function(model, contrast, SE) {
    if (class(model)[1] == "try-error") {
        return(NA)
    }
    if (nrow(model[[SE]]) == length(contrast)) {
        return(diag(t(contrast) %*% model[[SE]] %*% contrast))
    } else {
        return(NA)
    }
}

#' Perform a Wald test on GEE models
#'
#' @description Perform a Wald test on GEE models
#'
#' @param object A `SingleCellExperiment` object that already contains fitted
#' GEE models in the `rowData` slot of the object, as obtained with the GEEutils
#' function `bakerySCE`.
#'
#' @param contrast A `matrix` specifying the contrast, i.e. combinations of
#'   model parameters, of interest.
#'
#' @param SE A `character string` that indicates
#'
#' @return A `Dataframe` containing the requested model parameter estimates,
#'   (robust) standard errors, degrees of freedom, Wald test statistics and
#'   p-values.
#'
#' @author Jeroen Gilis
#'
#' @examples
#' library(scuttle)
#' set.seed(011235)
#' sce <- mockSCE(ncells = 400, ngenes = 10)
#' sce$patient_id <- factor(rep(paste0("patient", 1:8), each = ncol(sce) / 8))
#' sce$group_id <- factor(rep(paste0("group", 1:2), each = ncol(sce) / 2))
#'
#' metadata(sce)$formula <- ~ sce$group_id
#'
#' ## Fit models with bakerySCE
#' sce_fitted <- bakerySCE(
#'     object = sce,
#'     id = "patient_id",
#'     corstr = "independence",
#'     extraSandwich = "none",
#'     family = "poisson"
#' )
#'
#' ## Set up design and contrast
#' design <- model.matrix(~group_id, data = colData(sce))
#' L <- matrix(data = 0, ncol = 1, nrow = ncol(design)) # single contrast
#' rownames(L) <- colnames(design)
#' L["group_idgroup2", 1] <- 1
#' colnames(L) <- "group1_vs_group2"
#'
#' result <- geeWaldTest(object = sce_fitted, contrast = L[, 1], SE = "robust.variance")
#' ## result is a `data.frame` that contains the results of the Wald test for
#' ##  the specified contrasts and working under the robust standard error by
#' ## Liang and Zeger (1986).
#' @export
geeWaldTest <- function(object, contrast, SE) {
    models <- rowData(object)[["geefit"]]

    estimates <- vapply(models,
        FUN = .gee_getEstimates, FUN.VALUE = numeric(1),
        contrast = contrast
    )
    var <- vapply(models,
        FUN = .gee_varContrast, FUN.VALUE = numeric(1),
        contrast = contrast, SE = SE
    )

    se <- sqrt(var)
    W_stats <- estimates / se
    pvals <- 2 * pnorm(abs(W_stats), lower.tail = FALSE) # Z-stat
    data.frame(estimates, se, W_stats, pvals)
}
