% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterBySample.R
\name{filterBySample}
\alias{filterBySample}
\title{Filter genes by expression per sample}
\usage{
filterBySample(object, sample_id, size = 2, use_assay = 1L)
}
\arguments{
\item{object}{A \linkS4class{SingleCellExperiment} object}

\item{sample_id}{A \code{character} which identifies the samples or experimental
units. Must be one of the column names in the \code{colData} of the provided
\linkS4class{SingleCellExperiment} object}

\item{size}{An \code{integer} defining the minimum number of samples that must
express a gene for the gene to pass the filtering. Defaults to 2.}

\item{use_assay}{A \code{character} or \code{integer} specifying which \code{assay} from the
\linkS4class{SingleCellExperiment} to use. The default is to use the first
assay in \code{assayNames(object)}.}
}
\value{
An updated \linkS4class{SingleCellExperiment} object, with all genes
that did not pass the filtering removed.
}
\description{
Function for gene-level filtering. Removes genes that are expressed in less
than a certain number of samples (default = 2). Note: "samples" here refers
to the biological subjects from which the data was obtained.
}
\examples{
library(scuttle)
set.seed(011235)
sce <- mockSCE(ncells = 400, ngenes = 200)
sce$patient_id <- factor(rep(paste0("patient", 1:8), each = ncol(sce) / 8))

## Set some genes to 0 for all but 1 patient
counts(sce)[1:10, -which(sce$patient_id == "patient1")] <- 0

sce_filt <- filterBySample(sce, sample_id = "patient_id", size = 2)
sce_filt # filtered SCE object

}
\author{
Jeroen Gilis, Milan Malfait
}
