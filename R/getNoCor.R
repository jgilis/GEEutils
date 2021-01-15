#' Function to obtain a count matrix of uncorrelated features
#'
#' @description Retain only uncorrelated features from a count matrix
#'
#' @param object Either a `matrix` or a `SingleCellExperiment` instance. If it is a `SingleCellExperiment` instance, the first assay will be taken
#'  as the reference assay from which uncorrelated features are selected.
#'
#' @param nFeatures The number of correlated features to be returned. Note that if this number is larger than the number of uncorrelated features
#'  in the data, a `matrix` or `SingleCellExperiment` with as many uncorrelated features as possible will be returned. Default is 500 features.
#'
#' @param method A `charachter vector` indicating which correlation method should be used. One of "spearman", "pearson" or "kendall".
#'  Default is "spearman".
#'
#' @param cutoff A `numeric value` between 0 and 1, indicating the maximum allowed correlation between the returned features. Defaults to 0.25.
#'
#' @return A new `matrix` or `SingleCellExperiment` instance
#'
#' @rdname getNoCor
#'
#' @author Jeroen Gilis
#'
#' @export

getNoCor <- function(object, nFeatures = 500, method = "spearman", cutoff = 0.25){

  if(class(object)[1] == "matrix"){
      counts <- as.data.frame(t(object)) # transpose necessary for cor function
  } else if(class(object)[1] == "SingleCellExperiment"){
      counts <- as.data.frame(t(as.matrix(sce@assays@data[[1]]))) # transpose necessary for cor function
  } else{
    message("The provided object is not a matrix nor a SingleCellExperiment")
    stop()
  }

  all_corr <- cor(counts,method = method) # compute all pairwise correlations, so gene x gene matrix. If too big, one code sample 1 by 1 as Zimmerman.
  selected_genes <- c()

  for (i in seq_len(nFeatures)) {
    current_gene <- sample(colnames(all_corr),1)
    selected_genes <- c(selected_genes,current_gene)
    remove <- which(abs(all_corr[current_gene,]) > cutoff)
    all_corr <- all_corr[-remove,-remove]
    if (ncol(all_corr) < 10) break # if too few uncorrelated features left --> break
  }
  return(object[selected_genes,])
}
