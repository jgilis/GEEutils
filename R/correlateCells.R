# TODO: add option to specify number pairs to consider (sample from all possible pairs)
# TODO: add `type = "both"` to calculate both types in one go (return as list?)
#' @importFrom BiocParallel SerialParam
.correlate_cells <- function(x, grouping = NULL,
                             type = c("within", "between"),
                             BPPARAM = SerialParam()) {
    type <- match.arg(type)
    n_cells <- ncol(x)

    if (is.null(grouping)) {
        if (type == "between") {
            stop(
                "`type = 'between'` requires `grouping` to be specified",
                call. = FALSE
            )
        }
        cell_idx <- list(seq_len(n_cells))
    } else {
        if (length(grouping) != n_cells) {
            stop(
                "`grouping` should have the same length as the number of cells",
                call. = FALSE
            )
        }
        cell_idx <- split(seq_len(n_cells), grouping, drop = TRUE)
    }

    pairings <- .get_cell_pairs(cell_idx, type = type)
    .calculate_cor(x, pairings, BPPARAM = BPPARAM)
}


#' @importFrom utils combn
.get_cell_pairs <- function(cell_idx, type) {
    if (type == "within") {
        out <- lapply(cell_idx, combn, m = 2L)
    } else if (type == "between") {
        id_grid <- expand.grid(cell_idx, KEEP.OUT.ATTRS = FALSE)
        out <- apply(as.matrix(id_grid), 1, combn, m = 2L, simplify = FALSE)
    } else {
        stop("Invalid `type` argument.", call. = FALSE)
    }
    ## Output: 2-column matrix of cell ID pairs
    out <- t(do.call("cbind", out))
    unique(out) # remove any redundancies
}


.calculate_cor <- function(x, pairings, BPPARAM = BPPARAM) {
    ## Hack scran::correlatePairs() to correlate cells instead of genes
    out <- scran::correlatePairs(t(x), pairings = pairings, BPPARAM = BPPARAM)
    colnames(out)[c(1, 2)] <- c("cell1", "cell2")
    out
}


## S4 stuff
setGeneric("correlateCells", function(x, ...) standardGeneric("correlateCells"))

setMethod("correlateCells", "ANY", .correlate_cells)

setMethod("correlateCells", "SummarizedExperiment",
    function(x, grouping = NULL, ..., use_assay = "logcounts") {
        if (is.character(grouping) && length(grouping) == 1) {
            grouping <- colData(x)[[grouping]]
        }
        .correlate_cells(assay(x, use_assay), grouping = grouping, ...)
    }
)
