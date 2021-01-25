# Compute pairwise inter-individual correlation
#' @importFrom utils combn
.interpairwiseCor <- function(X) {
    pairs <- combn(colnames(as.data.frame(X)), 2, simplify = FALSE)
    ngenes_vector <- vector(mode = "numeric", length = length(pairs)) # not used?
    corr_vector <- vector(mode = "numeric", length = length(pairs))
    for (i in 1:length(pairs)) {
        cell1 <- as.numeric(X[, pairs[[i]][1]])
        cell2 <- as.numeric(X[, pairs[[i]][2]])
        corr_vector[i] <- round(cor(cell1, cell2, method = "spearman"), 4)
    }
    corr_vector
}

# Extract a single cell from an expression matrix
.pullonecell <- function(x) {
    x[sample(1:nrow(x), 1), ]
}

# Apply .pullonecell on an expression matrix
.onecell <- function(x, data) {
    # FIXME: replace `sapply` with `vapply`
    # cfr. https://bioconductor.org/developers/package-guidelines/#rcode
    sapply(data, .pullonecell)
}

# Compute estimates for target contrast
#' @importFrom utils combn
.intrapairwiseCor <- function(X) {
    pairs <- combn(rownames(as.data.frame(X)), 2, simplify = FALSE)
    pairs <- sample(pairs, 60) # comment out for all
    corr_vector <- vector(mode = "numeric", length = length(pairs))
    for (i in 1:length(pairs)) {
        cell1 <- as.numeric(X[pairs[[i]][1], ])
        cell2 <- as.numeric(X[pairs[[i]][2], ])
        corr_vector[i] <- round(cor(cell1, cell2, method = "spearman"), 4)
    }
    corr_vector
}


#' Visualize inter- and intra-individual correlations by means of boxplots
#'
#' @description Visualize inter- and intra-individual correlations by means of
#'   boxplots
#'
#' @param object A `SingleCellExperiment` object.
#'
#' @param clustering A `character` string indicating from which individual each
#'   cells originates
#'
#' @return A `ggplot` object.
#'
#' @rdname correlationsBoxplot
#'
#' @author Jeroen Gilis
#'
#' @import ggplot2
#'
#' @export
correlationsBoxplot <- function(object, clustering) {
    message("Function still under construction")

    # wrangle data
    # FIXME: use of `@` heavily discouraged; use appropriate accessor instead
    data <- as.data.frame(t(as.matrix(object@assays@data@listData$counts)))
    data <- split(data, object[[clustering]])

    # inter-individual correlations
    n_rep <- 10
    # TODO: I think this can be further optimizd by combining the 2 loops
    repmeans <- lapply(seq_len(n_rep), .onecell, data = data)
    repmeans_inter_corr <- lapply(repmeans, .interpairwiseCor)

    inter_corr <- gdata::unmatrix(do.call("rbind", repmeans_inter_corr))
    inter_corr <- cbind(
        rep("Inter_Alpha_Correlation", length(inter_corr)),
        inter_corr
    )
    inter_corr <- as.data.frame(
        cbind(rep("Inter", nrow(inter_corr)), inter_corr)
    )
    colnames(inter_corr) <- c("Type", "CellType", "Correlation")
    inter_corr$Correlation <- as.numeric(as.character(inter_corr$Correlation))

    # intra-individual correlations
    intra_corr <- as.numeric(unlist(sapply(data, .intrapairwiseCor)))
    intra_corr <- cbind(
        rep("Intra_Alpha_Correlation",
        length(intra_corr)),
        intra_corr
    )
    intra_corr <- as.data.frame(cbind(
        rep("Intra",
        nrow(intra_corr)),
        intra_corr
    ))
    colnames(intra_corr) <- c("Type", "CellType", "Correlation")
    intra_corr$Correlation <- as.numeric(as.character(intra_corr$Correlation))

    allcorr <- rbind(inter_corr, intra_corr)

    # boxplot
    ggplot(allcorr, aes(x = CellType, y = Correlation, fill = Type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(pch = 21, position = position_jitterdodge()) +
        theme_classic() +
        # scale_y_continuous(limits=c(0.2,0.6),breaks=(seq(-10,10,1)/10)) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
}
