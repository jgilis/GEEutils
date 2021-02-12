#' Compute robust standard error estimates using multiple sandwich strategies
#' for SingleCellExperiment input
#'
#' @description Compute robust standard error estimates using multiple sandwich
#'   strategies for SingleCellExperiment input
#'
#' @param object A `SingleCellExperiment` object
#'
#' @param id A `vector` which identifies the clusters. The length of id should
#'   be the same as the number of observations. Data are assumed to be sorted so
#'   that observations on a cluster are contiguous rows for all entities in the
#'   formula.
#'
#' @param corstr A `character string` specifying the correlation structure.
#'   Default is "exchangeable". The following are permitted: "independence",
#'   "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and
#'   "unstructured".
#'   
#' @param extraSandwich A `character vector` of sandwich estimation procedures
#'   that must be used to compute robust standard errors, additional to the
#'   model standard error and the robust estimator from Liang and Zeger. The
#'   following are permitted: "none", "DF", "KC" and "Pan". If not specified or 
#'   set to "none", only the naive model standard error and the robust estimator 
#'   from Liang and Zeger will be provided.
#'
#' @param family A `character string` indicating the family for defining link
#'   and variance functions. Currently, only "poisson" is supported.
#'
#' @return An updated `SingleCellExperiment` object. The fitted GEE models for
#' each gene can be retrieved from the `rowData` slot.
#'
#' @rdname bakerySCE
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
#' ## Run bakery with KC and Pan extra sandwiches
#' out <- bakery(
#'     formula = gene ~ group_id,
#'     id = "patient_id",
#'     data = data,
#'     family = "poisson",
#'     corstr = "exchangeable",
#'     extraSandwich = c("KC", "Pan")
#' )
#'
#' ## Liang-Zeger variance estimates
#' diag(out$robust.variance)
#'
#' ## KC variance estimates
#' diag(out$KC.variance)
#'
#' ## Pan variance estimates
#' diag(out$Pan.variance)
#'
#' @importFrom stats model.matrix
#' @importFrom MultiAssayExperiment metadata
#' @importFrom SummarizedExperiment assays colData rowData rowData<-
#'
#' @export
bakerySCE <- function(object, 
                      id, 
                      corstr, 
                      extraSandwich = "none", 
                      family = "poisson") {
  
  counts <- assays(object)[["counts"]] #hard-coded to take assay "counts"
  data <- as.data.frame(colData(object))
  design <- model.matrix(metadata(object)$formula, data = data)
  
  geefit <- lapply(1:nrow(counts), function(i) {
    
    #TODO: could wrap an updated eval_fork around bakery
    geefit_i <- suppressMessages(bakery(
      formula = counts[i,] ~ -1+design,
      id = id,
      data = data,
      family = family,
      corstr = corstr,
      extraSandwich = extraSandwich
    ))
  })
  names(geefit) <- rownames(counts)
  
  rowData(object)[["geefit"]] <- geefit # store in rowData
  return(object)
}