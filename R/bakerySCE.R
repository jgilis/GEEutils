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
#' set.seed(011235)
#' sce <- scuttle::mockSCE(ncells = 400, ngenes = 10)
#' sce$patient_id <- factor(rep(paste0("patient", 1:8), each = ncol(sce) / 8))
#' sce$group_id <- factor(rep(paste0("group", 1:2), each = ncol(sce) / 2))
#'
#' MultiAssayExperiment::metadata(sce)$formula <- ~sce$group_id
#'
#' sce_fitted <- bakerySCE(object = sce,
#'                         id = "patient_id",
#'                         corstr = "independence",
#'                         extraSandwich = "none",
#'                         family = "poisson") 
#'
#' ## Model fit for gene "Gene_0001"
#' SummarizedExperiment::rowData(sce_fitted)[["geefit"]]$Gene_0001
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