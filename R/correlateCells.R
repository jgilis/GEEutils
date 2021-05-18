# TODO: add `BPPARAM` argument to pass on to `scran::correlatePairs()`
.correlate_cells <- function(x, grouping = NULL, type = "within") {
    type <- match.arg(type)

    n_cells <- ncol(x)

    if (is.null(grouping)) {
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
    .calculate_cor(x, pairings)
}


#' @importFrom utils combn
.get_cell_pairs <- function(cell_idx, type) {
    if (type == "within") {
        out <- lapply(cell_idx, combn, m = 2L)
    }
    out <- do.call("cbind", out)
    ## Output: 2-column matrix of cell ID pairs
    t(out)
}


.calculate_cor <- function(x, pairings) {
    ## Hack scran::correlatePairs() to correlate cells instead of genes
    out <- scran::correlatePairs(t(x), pairings = pairings)
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
