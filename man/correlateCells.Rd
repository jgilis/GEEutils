% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlateCells.R
\name{correlateCells}
\alias{correlateCells}
\alias{correlateCells,ANY-method}
\alias{correlateCells,SummarizedExperiment-method}
\alias{plotCorrelations}
\alias{plotCorrelations,DataFrame-method}
\alias{plotCorrelations,SummarizedExperiment-method}
\title{Compare within- and between-group correlations between cells}
\usage{
correlateCells(x, ...)

\S4method{correlateCells}{ANY}(x, grouping, n_pairs = NULL, BPPARAM = SerialParam())

\S4method{correlateCells}{SummarizedExperiment}(x, grouping, ..., use_assay = "logcounts")

plotCorrelations(x, ...)

\S4method{plotCorrelations}{DataFrame}(x, point_size = 0.5)

\S4method{plotCorrelations}{SummarizedExperiment}(x, ..., point_size = 0.5)
}
\arguments{
\item{x}{For \code{correlateCells}, a numeric \code{matrix} of counts where genes are
rows and cells are columns. Alternatively, a
\linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment}
object.

For \code{plotCorrelations}, a \code{DataFrame} containing the
output from \code{correlateCells}. Alternatively, a
\linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment}
object, in which case \code{correlateCells} will be run on the fly.}

\item{...}{For the generics, further arguments passed to specific methods.

For the \code{correlateCells} \emph{SummarizedExperiment} method, further arguments
to pass to the \emph{ANY} method.

For the \code{plotCorrelations} \emph{SummarizedExperiment} method, further arguments
to pass to \code{correlateCells}.}

\item{grouping}{A character or numeric vector of length \code{ncol(x)} specifying
the group membership of the cells. If \code{x} is a \emph{SummarizedExperiment} or
\emph{SingleCellExperiment}, can also pass the name of a column in the
\code{colData}.}

\item{n_pairs}{Optional integer to specify the number of cell pairs used to
calculate correlations, If \code{NULL} (the default), will use all possible
pairs. Note that this can become intractable fairly quickly for larger data
sets. Setting \code{n_pairs} will result in randomly sub-sampling of the
within-group and between-group pairs such that they are both roughly equal
to \code{n_pairs} (so the total number pairs will be roughly \code{2 * n_pairs}).
Note that deviations are likely due to the sampling of redundant pairs,
which are then removed.}

\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying whether
the correlation calculations should be parallelized.}

\item{use_assay}{If \code{x} is a \emph{SummarizedExperiment}, a character string
specifying which assay to use from \code{x}. The default is to use
\code{"logcounts"}. It's recommended to use log-normalized counts for the
correlations.}

\item{point_size}{Numeric scalar, specifying the size of the points in the
plot.}
}
\value{
For \code{correlateCells}, a \code{DataFrame} with one row for each pair of
cells with the following columns:
\itemize{
\item \code{cell1}, \code{cell2}: character strings specifying the names of the cells in
the pair
\item \code{rho}: estimate of Spearman's correlation (\eqn{\rho}).
\item \code{p.value}, \code{FDR}: the p-value and its FDR-corrected equivalent for the null
hypothesis of no correlation between the cells
\item \code{type}: character string denoting the type of correlation, either
\strong{within-group} or \strong{between-group}.
}

For \code{plotCorrelations}, a \code{ggplot} object.
}
\description{
Calculate and plot correlations between cells from the same group and between
cells from different groups, based on Spearman's rank correlation. The groups
here are typically subjects or samples, to diagnose whether correlations tend
to be higher for cells within the same biological unit.
}
\details{
This is mainly meant as a diagnostic tool to investigate whether there tends
to be a higher correlation between cells from the same biological unit (e.g.
subject, sample, patient, ...). If this is the case, this should be taken
into account when performing differential expression (DE) analyses. For
example by using (clustered) robust variance estimators, as implemented in
\code{\link[=glmSandwichTest]{glmSandwichTest()}}.
}
\examples{
library(scuttle)
n_groups <- 3L
n_cells <- 50L

sce <- mockSCE(ncells = n_cells, ngenes = 100)
sce$subject <- sample(LETTERS[1:n_groups], n_cells, replace = TRUE)
sce <- logNormCounts(sce)

## Calculate within- and between-subject correlations
corr <- correlateCells(sce, grouping = "subject", n_pairs = 100)
head(corr)

## Plot previously calculated correlations
plotCorrelations(corr)

## ...or directly run plotCorrelations() on the sce object
plotCorrelations(sce, grouping = "subject", n_pairs = 100)

}
\seealso{
\code{\link[scran:correlatePairs]{scran::correlatePairs()}} which is used under the hood, but on the cells
instead of the genes.

\code{\link[=glmSandwichTest]{glmSandwichTest()}} for a DE method that can take into
account the correlation structure.
}
\author{
Milan Malfait
}
