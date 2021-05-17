.correlate_cells <- function(x, grouping = NULL, type = "within") {
    type <- match.arg(type)
    cor_fun <- switch(type,
        within = .within_cor
    )

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

    cor_fun(x, cell_idx)
}


.get_cell_pairs <- function(cell_idx) {
    out <- combn(cell_idx, 2L)
    list(cell1 = out[1, ], cell2 = out[2, ])
}


## Within-group cell-wise correlations
.within_cor <- function(x, cell_idx) {
    rhos <- vector("list", length = length(cell_idx))

    for (i in seq_along(rhos)) {
        current_ids <- cell_idx[[i]]
        pairs <- .get_cell_pairs(current_ids)
        x1 <- x[, pairs$cell1, drop = FALSE]
        x2 <- x[, pairs$cell2, drop = FALSE]

        rhos[[i]] <- diag(cor(x1, x2, method = "spearman"))
    }
    unlist(rhos)
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
