# Wrapper to prettify reduced dimension plots (from the Muscat R package)
.plot_dr <- function(object, type, clustering) {
      scater::plotReducedDim(object, dimred = type, colour_by = clustering) +
          guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
          theme_minimal() +
          theme(aspect.ratio = 1)
  }

# TODO: do we really need this as a package function?
# Seems a bit out of line with rest of pkg; unless we would add more plotting functionality

#' Wrapper function to produce a pretty dimensionality reduction plot
#'
#' @description Wrapper function to produce a pretty dimensionality reduction
#'   plot
#'
#' @param object A `SingleCellExperiment` instance. Must have 2D dimensionality
#'   reduced coordinates precomputed, which can be checked with
#'   `reducedDimNames(object)`.
#'
#' @param clustering A `character vector` indicating the categorical variable
#'   that is used to color the dimensionality reduction plot. Must be one of
#'   `colnames(colData(object))`.
#'
#' @param type A `character vector` indicating which dimensionality reduction
#'   method should be used for visualization. Must be one of
#'   `reducedDimNames(object)`. Defaults to "TSNE".
#'
#'
#' @param downsample Should there be a downsampling of the number of cells per
#'   cluster in the plot? Either `logical` FALSE if no downsampling is desired
#'   or a `numeric` indicating the maximum number of cells per cluster.
#'
#' @return Prints the dimensionality reduced plot directly to console.
#'
#' @rdname runDimRed
#'
#' @author Jeroen Gilis
#'
#' @importFrom scater plotReducedDim
#' @import ggplot2
#'
#' @export
# TODO: why not return as ggplot object? Would allow further manipulation and
#       exporting of plot
runDimRed <- function(object, clustering, type = "TSNE", downsample = 500) {
    if (downsample) {
        # downsample to max. "downsample" cells per cluster
        cs_by_k <- split(colnames(object), clustering)
        # FIXME: replace `sapply` with `vapply`
        # cfr. https://bioconductor.org/developers/package-guidelines/#rcode
        csDown <- unlist(sapply(cs_by_k, function(u) {
              sample(u, min(length(u), downsample))
          }))

        # plot t-SNE
        print(.plot_dr(object[, csDown], type, clustering))
    } else {
        # plot t-SNE
        print(.plot_dr(object, type, clustering))
    }
}
