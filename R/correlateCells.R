# TODO: add plotting function
#' @importFrom BiocParallel SerialParam
.correlate_cells <- function(x, grouping,
                             n_pairs = NULL,
                             BPPARAM = SerialParam()) {
    n_cells <- ncol(x)

    if (is.null(grouping)) {
        ## Assume all cells from single group
        grouping <- 1
    }
    cell_idx <- split(seq_len(n_cells), grouping, drop = TRUE)

    pairings <- .get_cell_pairs(cell_idx, n_pairs = n_pairs)
    all_pairs <- rbind(pairings$within, pairings$between)

    out <- .calculate_cor(x, pairings = all_pairs, BPPARAM = BPPARAM)
    out[["type"]] <- c(
        rep("within", NROW(pairings$within)),
        rep("between", NROW(pairings$between))
    )
    out
}


#' @importFrom utils combn
#' @importFrom Matrix t
.get_cell_pairs <- function(cell_idx, n_pairs) {
    within <- .within_pairs(cell_idx, n_pairs = n_pairs)
    if (length(cell_idx) > 1) {
        between <- .between_pairs(cell_idx, n_pairs = n_pairs)
    } else {
        between <- NULL
    }
    list(within = within, between = between)
}

.within_pairs <- function(cell_idx, n_pairs) {
    ## Omit subjects with just 1 cell; combn() uses seq_len(n) for integer n
    cell_idx[lengths(cell_idx) == 1] <- NULL
    out <- lapply(cell_idx, combn, m = 2L)

    ## Output: 2-column matrix of cell ID pairs
    out <- t(do.call("cbind", out))
    out <- unique(out) # remove any redundancies

    if (!is.null(n_pairs)) {
        if (n_pairs > nrow(out)) {
            warning(
                "`n_pairs` greater than max. possible within-group pairs.",
                "\n  Will use all pairs (n = ", nrow(out), ")."
            )
        } else {
            out <- out[sample(nrow(out), size = n_pairs), ]
        }
    }
    out
}

.between_pairs <- function(cell_idx, n_pairs) {
    if (is.null(n_pairs)) {
        ## Use all possible between-group pairs
        id_grid <- expand.grid(cell_idx, KEEP.OUT.ATTRS = FALSE)
        out <- apply(as.matrix(id_grid), 1, combn, m = 2L, simplify = FALSE)
    } else {
        ## Use sampling strategy to avoid computing all possible pairs
        groups <- names(cell_idx)
        group_pairs <- t(combn(groups, 2))

        ## Sample a number of cells per group-pair such that we get n_pairs total
        cells_per_pair <- round(n_pairs / nrow(group_pairs))

        out <- apply(group_pairs, 1, function(g) {
            g1 <- g[[1]]
            g2 <- g[[2]]
            rbind(
                sample(cell_idx[[g1]], size = cells_per_pair, replace = TRUE),
                sample(cell_idx[[g2]], size = cells_per_pair, replace = TRUE)
            )
        }, simplify = FALSE)
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
