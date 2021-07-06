#' Compare within- and between-group correlations between cells
#'
#' Calculate and plot correlations between cells from the same group and between
#' cells from different groups, based on Spearman's rank correlation. The groups
#' here are typically subjects or samples, to diagnose whether correlations tend
#' to be higher for cells within the same biological unit.
#'
#' @param x For `correlateCells`, a numeric `matrix` of counts where genes are
#'   rows and cells are columns. Alternatively, a
#'   \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment}
#'   object.
#'
#'   For `plotCorrelations`, a `DataFrame` containing the
#'   output from `correlateCells`. Alternatively, a
#'   \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment}
#'   object, in which case `correlateCells` will be run on the fly.
#'
#' @param grouping A character or numeric vector of length `ncol(x)` specifying
#'   the group membership of the cells. If `x` is a *SummarizedExperiment* or
#'   *SingleCellExperiment*, can also pass the name of a column in the
#'   `colData`.
#' @param n_pairs Optional integer to specify the number of cell pairs used to
#'   calculate correlations, If `NULL` (the default), will use all possible
#'   pairs. Note that this can become intractable fairly quickly for larger data
#'   sets. Setting `n_pairs` will result in randomly sub-sampling of the
#'   within-group and between-group pairs such that they are both roughly equal
#'   to `n_pairs` (so the total number pairs will be roughly `2 * n_pairs`).
#'   Note that deviations are likely due to the sampling of redundant pairs,
#'   which are then removed.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   the correlation calculations should be parallelized.
#' @param ... For the generics, further arguments passed to specific methods.
#'
#'   For the `correlateCells` *SummarizedExperiment* method, further arguments
#'   to pass to the *ANY* method.
#'
#'   For the `plotCorrelations` *SummarizedExperiment* method, further arguments
#'   to pass to `correlateCells`.
#'
#' @param use_assay If `x` is a *SummarizedExperiment*, a character string
#'   specifying which assay to use from `x`. The default is to use
#'   `"logcounts"`. It's recommended to use log-normalized counts for the
#'   correlations.
#' @param point_size Numeric scalar, specifying the size of the points in the
#'   plot.
#'
#' @details
#' This is mainly meant as a diagnostic tool to investigate whether there tends
#' to be a higher correlation between cells from the same biological unit (e.g.
#' subject, sample, patient, ...). If this is the case, this should be taken
#' into account when performing differential expression (DE) analyses. For
#' example by using (clustered) robust variance estimators, as implemented in
#' [glmSandwichTest()].
#'
#' @return
#' For `correlateCells`, a `DataFrame` with one row for each pair of
#' cells with the following columns:
#'
#' * `cell1`, `cell2`: character strings specifying the names of the cells in
#'   the pair
#' * `rho`: estimate of Spearman's correlation (\eqn{\rho}).
#' * `p.value`, `FDR`: the p-value and its FDR-corrected equivalent for the null
#'   hypothesis of no correlation between the cells
#' * `type`: character string denoting the type of correlation, either
#'   __within-group__ or __between-group__.
#'
#' For `plotCorrelations`, a `ggplot` object.
#'
#' @author Milan Malfait
#' @name correlateCells
#'
#' @seealso
#' [scran::correlatePairs()] which is used under the hood, but on the cells
#' instead of the genes.
#'
#' [glmSandwichTest()] for a DE method that can take into
#' account the correlation structure.
#'
#' @examples
#' library(scuttle)
#' n_groups <- 3L
#' n_cells <- 50L
#'
#' sce <- mockSCE(ncells = n_cells, ngenes = 100)
#' sce$subject <- sample(LETTERS[1:n_groups], n_cells, replace = TRUE)
#' sce <- logNormCounts(sce)
#'
#' ## Calculate within- and between-subject correlations
#' corr <- correlateCells(sce, grouping = "subject", n_pairs = 100)
#' head(corr)
#'
#' ## Plot previously calculated correlations
#' plotCorrelations(corr)
#'
#' ## ...or directly run plotCorrelations() on the sce object
#' plotCorrelations(sce, grouping = "subject", n_pairs = 100)
#'
NULL

#' @importFrom BiocParallel SerialParam
.correlate_cells <- function(x, grouping,
                             n_pairs = NULL,
                             BPPARAM = SerialParam()) {
    n_cells <- ncol(x)

    if (is.null(grouping)) {
        ## Assume all cells from single group
        grouping <- 1
    } else {
        stopifnot(length(grouping) == ncol(x))
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
            out <- out[sample(nrow(out), size = n_pairs), , drop = FALSE]
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

        if (n_pairs < nrow(group_pairs)) {
            ## So that cells_per_pair does not become 0 in next step
            sub_pairs <- sample(nrow(group_pairs), n_pairs)
            group_pairs <- group_pairs[sub_pairs, , drop = FALSE]
        }

        ## Sample cells per group-pair such that we get n_pairs total
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


#' @rdname correlateCells
#' @export
setGeneric("correlateCells", function(x, ...) standardGeneric("correlateCells"))

#' @rdname correlateCells
#' @export
setMethod("correlateCells", "ANY", .correlate_cells)

#' @rdname correlateCells
#' @export
setMethod("correlateCells", "SummarizedExperiment",
    function(x, grouping, ..., use_assay = "logcounts") {
        if (is.character(grouping) && length(grouping) == 1) {
            grouping <- colData(x)[[grouping]]
        }
        .correlate_cells(assay(x, use_assay), grouping = grouping, ...)
    }
)


# Plotting ----------------------------------------------------------------

#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point position_jitter
.plot_correlations <- function(x, point_size = 0.5) {
    df <- as.data.frame(x)
    type <- rho <- NULL
    ggplot(df, aes(x = type, y = rho)) +
        geom_boxplot(aes(fill = type), outlier.size = -1, alpha = 0.5) +
        geom_point(size = point_size, position = position_jitter(width = 0.2))
}


#' @rdname correlateCells
#' @aliases plotCorrelations
#' @export
setGeneric("plotCorrelations",
    function(x, ...) standardGeneric("plotCorrelations")
)

## If DataFrame is given, assume this is the output from correlateCells()
#' @rdname correlateCells
#' @aliases plotCorrelations
#' @export
setMethod("plotCorrelations", "DataFrame", .plot_correlations)

## If SummarizedExperiment given, first run correlateCells() before plotting
#' @rdname correlateCells
#' @aliases plotCorrelations
#' @export
setMethod("plotCorrelations", "SummarizedExperiment",
          function(x, ..., point_size = 0.5) {
              corr <- correlateCells(x, ...)
              plotCorrelations(corr, point_size = point_size)
          }
)
