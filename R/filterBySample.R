#' Function for gene-level filtering.
#'
#' @description Function for gene-level filtering. Removes genes that are
#' expressed in less than a certain number of samples (default = 2) in any of
#' the experimental groups of interest. To be used in combination with
#' `filterByExpr` from the edgeR package.
#'
#' @param object A `SingleCellExperiment` object
#'
#' @param sample_id A `vector` which identifies the samples or experimental
#' units. Must be one of the column names in the `colData` of the provided
#' `SingleCellExperiment` object
#'
#' @param size A `number` defining the minimum amount of samples that must
#' express a gene for the gene to pass the filtering. Defaults to 2.
#'
#' @param use_assay A `character` or `integer` specifying which `assay` from the
#' `SingleCellExperiment` to use. Default: `1L`, which will use the first assay
#' in `assayNames(object)`.
#'
#' @return An updated `SingleCellExperiment` object, with all genes that did
#' not pass the filtering removed.
#'
#' @author Jeroen Gilis
#'
#' @examples
#' library(scuttle)
#' set.seed(011235)
#' sce <- mockSCE(ncells = 400, ngenes = 200)
#' sce$patient_id <- factor(rep(paste0("patient", 1:8), each = ncol(sce) / 8))
#'
#' ## Set some genes to 0 for all but 1 patient
#' counts(sce)[1:10, -which(sce$patient_id == "patient1")] <- 0
#'
#' sce_filt <- filterBySample(sce, sample_id = "patient_id", size = 2)
#' sce_filt # filtered SCE object
#'
#' @export
filterBySample <- function(object, sample_id,
                           size = 2, use_assay = 1L) {

    # Binarize assay -> upon aggregation, we count the number of cells
    # with any expression for each gene
    binary_assay <- assay(object, use_assay) > 0

    samples <- factor(object[[sample_id]])
    split_cells <- split(seq_len(ncol(object)), f = samples, drop = TRUE)

    pb <- .counts_per_sample(x = binary_assay, split_cells = split_cells)
    keep <- rowSums(pb > 0) > size
    # NOTE: the >0 could also be >count, with count an input of filterBySample
    object[keep, ]
}


#' @importFrom Matrix rowSums
.counts_per_sample <- function(x, split_cells) {
    out <- lapply(split_cells, function(i) {
        if (length(i) == 0) {
            return(numeric(nrow(x)))
        }
        rowSums(x[, i, drop = FALSE])
        # set aggregation function to Matrix::rowSums <-> muscat
    })
    do.call(cbind, out)
}
