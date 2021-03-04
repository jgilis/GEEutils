# Adapted from Muscat
# See: https://github.com/HelenaLC/muscat/, utils.R 
# ------------------------------------------------------------------------------
# split cells by group-sample
# ------------------------------------------------------------------------------
#   x:  a SingleCellExperiment
#   by: character vector specifying colData column(s) to split by
# ------------------------------------------------------------------------------
#' @importFrom SummarizedExperiment colData
#' @importFrom data.table data.table
#' @importFrom purrr map_depth
.split_cells_filter <- function (x, by)
{
  if (is(x, "SingleCellExperiment"))
    x <- colData(x)
  cd <- data.frame(x[by], check.names = FALSE)
  cd <- data.table(cd, cell = rownames(x))
  cd <- split(cd,
              by = by, 
              sorted = TRUE, 
              flatten = FALSE, 
              drop = TRUE) 
  # added drop=TRUE in split to remove empty group*sample combinations
  map_depth(cd, length(by), "cell")
}

# Adapted from Muscat
# See: https://github.com/HelenaLC/muscat/, utils-pbDS.R 
# ------------------------------------------------------------------------------
# Aggregate cells that originate from the same group-sample combination
# ------------------------------------------------------------------------------
#   x:  a SingleCellExperiment
#   cs: character vector specifying colData column(s) to split by
#   assay: defaults to first assay of the SingleCellExperiment
# ------------------------------------------------------------------------------
#' @importFrom Matrix rowSums
#' @importFrom purrr map_depth
#' @importFrom SummarizedExperiment assays
.pb_filter <- function (x, cs, assay) # function from muscat
{
  pb <- map_depth(cs, -1, function(i) {
    if (length(i) == 0) 
      return(numeric(nrow(x)))
    rowSums(assays(x)[[assay]][, i, drop = FALSE])
    # set aggregation function to Matrix::rowSums <-> muscat
  })
  map_depth(pb, -2, function(u) as.matrix(data.frame(u, 
                                                     row.names = rownames(x), 
                                                     check.names = FALSE)))
}

#' Function for gene-level filtering.
#' 
#' @description Function for gene-level filtering. Removes genes that are
#' expressed in less than a certain number of samples (default = 2) in any of 
#' the experimental groups of interest. To be used in combination with 
#' `filterByExpr` from the edgeR package.
#'
#' @param object A `SingleCellExperiment` object
#'
#' @param group_id A `vector` which identifies the experimental groups of 
#' interest. Must be one of the column names in the `colData` of the provided
#' `SingleCellExperiment` object
#'
#' @param sample_id A `vector` which identifies the samples or experimental 
#' units. Must be one of the column names in the `colData` of the provided
#' `SingleCellExperiment` object
#'
#' @param size A `number` defining the minimum amount of samples that must
#' express a gene for the gene to pass the filtering. Defaults to 2.
#'
#' @param assay A `character` specifying which `assay` from the
#' `SingleCellExperiment` to use. If `NULL` will use the first one in
#' `assayNames(object)`.
#'
#' @return An updated `SingleCellExperiment` object, with all genes that did
#' not pass the filtering removed.
#'
#' @rdname filterBySample
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
#' metadata(sce)$formula <- ~sce$group_id
#'
#' sce_filt <- filterBySample(object = sce,
#'                         group_id = "group_id",
#'                         sample_id = "patient_id",
#'                         size = 2,
#'                         assay = NULL)
#' sce_filt # filtered SCE object                  
#'
#' @importFrom SummarizedExperiment assays<- assayNames
#'
#' @export

filterBySample <- function(object, 
                           group_id, 
                           sample_id, 
                           size = 2,
                           assay = NULL){
  
  if (is.null(assay)) {
    assay <- assayNames(object)[[1]]
  }
  
  object_binary <- object # make copy
  assays(object_binary)[[assay]][assays(object_binary)[[assay]] > 0] <- 1
  # set counts > 0 to 1 -> upon aggregation, we count the number of cells
  # with any expression for each gene
  
  cs <- .split_cells_filter(x = object_binary, 
                     by = c(group_id, sample_id))
  
  pb <- .pb_filter(x = object_binary, 
                  cs = cs, 
                  assay = assay)
  
  remove <- c()
  for (i in seq_along(pb)) { # loop over group_id
    # remove genes if less than "size" samples have expression in 
    # a certain group
    remove <- c(names(which(rowSums(pb[[i]] > 0) <= size)), remove)
    # NOTE: the >0 could also be >count, with count an input of filterBySample
  }
  
  return(object[-which(rownames(object) %in% unique(remove)),])
}

