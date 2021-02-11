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
#' # TODO update example to SCE
#' ## Simulate single gene across 10 patients in 2 groups
#' n_cells <- 500
#' gene <- rpois(n_cells, 3)
#' group_id <- factor(rep(c("control", "treat"), each = n_cells / 2))
#' patient_id <- factor(rep(paste0("patient", 1:10), each = n_cells / 10))
#' data <- data.frame(gene, group_id, patient_id)
#'
#' @export
gee.Waldtest <- function(object, contrast, SE){
  # FIXME: replace `sapply` with `vapply`
  # cfr. https://bioconductor.org/developers/package-guidelines/#rcode
  # TODO: convert 'models' to list internally? Makes using a single model
  # bit more intuitive
  
  models <- rowData(object)[["fitDTUModels"]]
  
  estimates <- sapply(models, .gee_getEstimates, contrast = contrast)
  var <- sapply(models, .gee_varContrast, contrast = contrast, SE = SE)
  
  W_stats <- estimates/sqrt(var)
  pvals <- 2 * pnorm(abs(W_stats), lower.tail = FALSE) # Z-stat
  data.frame(estimates, var, W_stats, pvals)
}

