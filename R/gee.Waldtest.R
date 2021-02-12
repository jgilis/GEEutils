.gee_getCoef <- function(model){
  if(class(model)[1] == "try-error"){
    return(NULL)
  }
  return(model$coefficients)
}

.gee_getEstimates <- function(model, contrast){
  coef <- .gee_getCoef(model)
  if (is.null(coef)) {
    coef <- rep(NA, times = length(contrast))
  }
  contrast %*% coef
}

.gee_varContrast <- function(model,contrast,SE){
  if(class(model)[1] == "try-error"){
    return(NA)
  }
  if (nrow(model[[SE]]) == length(contrast)){
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
#' @rdname glm.Waldtest
#'
#' @author Jeroen Gilis
#'
#' @examples
#' set.seed(011235)
#' sce <- scuttle::mockSCE(ncells = 400, ngenes = 10)
#' sce$patient_id <- factor(rep(paste0("patient", 1:8), each = ncol(sce) / 8))
#' sce$group_id <- factor(rep(paste0("group", 1:2), each = ncol(sce) / 2))
#'
#' MultiAssayExperiment::metadata(sce)$formula <- ~sce$group_id
#'
#' ## Fit models with bakerySCE
#' sce_fitted <- bakerySCE(object = sce,
#'                         id = "patient_id",
#'                         corstr = "independence",
#'                         extraSandwich = "none",
#'                         family = "poisson") 
#'
#' ## Set up design and contrast
#' design <- stats::model.matrix(~group_id, data = SummarizedExperiment::colData(sce))
#' L <- matrix(data = 0, ncol = 1, nrow = ncol(design)) # single contrast
#' rownames(L) <- colnames(design)
#' L["group_idgroup2",1] <- 1
#' colnames(L) <- "group1_vs_group2"
#' 
#' result <- gee.Waldtest(object=sce_fitted, contrast=L[,1], SE="robust.variance")
#' ## result is a `data.frame` that contains the results of the Wald test for
#' ##  the specified contrasts and working under the robust standard error by
#' ## Liang and Zeger (1986).
#'
#' @export
gee.Waldtest <- function(object, contrast, SE){
  # FIXME: replace `sapply` with `vapply`
  # cfr. https://bioconductor.org/developers/package-guidelines/#rcode
  # TODO: convert 'models' to list internally? Makes using a single model
  # bit more intuitive
  
  models <- rowData(object)[["geefit"]]
  
  estimates <- sapply(models, .gee_getEstimates, contrast = contrast)
  var <- sapply(models, .gee_varContrast, contrast = contrast, SE = SE)
  
  W_stats <- estimates/sqrt(var)
  pvals <- 2 * pnorm(abs(W_stats), lower.tail = FALSE) # Z-stat
  data.frame(estimates, var, W_stats, pvals)
}

