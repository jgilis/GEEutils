#' Simulate a synthetic scRNA-Seq dataset based on parameters obtained with getSimulationParameters
#'
#' @description Simulate a synthetic scRNA-Seq dataset based on parameters obtained with getSimulationParameters
#'
#' @param design A `data frame` containing 3 columns.
#'  Column `cluster` is a `factor` that specifies the clustering of the observations.
#'  Column `group` is a `factor` that specifies to which (treatment) group the observations belong.
#'  Column `ncells` is a `numeric`indicating the number of cells belonging to each of the clusters.
#'
#' @param SimulationParameters A `list` of simulation parameters obtained with the `GEEutils::getSimulationParameters`.
#'
#' @param data An optional `data frame` in which to interpret the variables occurring in the formula,
#'  along with the id and n variables.
#'
#' @param nGenes A `numeric` indicating hoe many genes should be simulated.
#'
#' @param distribution A `character` indicating from which baseline distribution the data should be sampled. Either "NB" or "poisson".
#'
#' @param fractionDGE A `numeric` between 0 and 1, indicating the fraction of differentially expressed genes in the simulated data.
#' Defaults to zero.
#' @param foldchange A `numeric` indicating the size of differential expression in terms of a fold change.
#' Defaults to 1 (not differentially expressed, mock dataset).
#'
#' @return A `matrix` of simulated gene expression data.
#'
#' @rdname simulator
#'
#' @author Jeroen Gilis
#'
#' @import stats
#'
#' @export
simulator <- function(design,
                      SimulationParameters,
                      nGenes, distribution,
                      fractionDGE = 0,
                      foldchange = 1) {
    if (!"cluster" %in% colnames(design)) {
        stop("design must contain column with name cluster")
    }
    if (!"group" %in% colnames(design)) {
        stop("design must contain column with name group")
    }
    if (!"ncells" %in% colnames(design)) {
        stop("design must contain column with name ncells")
    }

    # TODO: might be more efficient to keep this as a matrix
    simulatedData <- as.data.frame(
        matrix(data = NA, nrow = nGenes, ncol = sum(design$ncells))
    )

    params <- SimulationParameters
    for (i in seq_len(nGenes)) {
        # get grandmean and relationship with inter-individual standard deviations
        grandmean <- rgamma(
            n = 1, shape = params$grandmean.shape, rate = params$grandmean.rate
        )
        stddev.of.within.means <- params$grandmean_interSD.beta1 * grandmean

        # get foldchange direction
        if (i <= nGenes * fractionDGE) { # its a DE gene
            # Now all foldchanges will be "foldchange" or "-foldchange".
            # TODO: It would be better to sample foldchanges from a distribution with mean "foldchange"
            foldchange <- ifelse(
                rbinom(n = 1, size = 1, prob = 0.5) == 1,
                foldchange,
                1 / foldchange
            )
        } else { # its not a DE gene
            foldchange <- 1
        }

        # get zero probability
        prob.zero <- rgamma(
            n = 1, shape = params$drop.shape, rate = params$drop.rate
        ) # sampled from gamma
        prob.zero <- ifelse(
            prob.zero > 1,
            rgamma(n = 1, shape = params$drop.shape, rate = params$drop.rate),
            prob.zero
        )
        drop.sd <- params$drop.beta0 + params$drop.beta1 * prob.zero + params$drop.beta2 * (prob.zero**2)
        drop.sd <- ifelse(drop.sd < 0, 0, drop.sd)
        prob.zero <- rnorm(n = 1, mean = prob.zero, sd = drop.sd)
        prob.zero <- ifelse(prob.zero < 0, 0, prob.zero)
        prob.zero <- ifelse(prob.zero > 1, 1, prob.zero)
        prob.zero <- 1 - prob.zero

        for (cluster_id in levels(design$cluster)) {
            # for each design levels of group
            reference <- design$group[1]

            if (design[which(design$cluster == cluster_id), "group"] == reference) {
                temp.mean <- grandmean + rnorm(n = 1, mean = 0, sd = stddev.of.within.means)
            } else {
                temp.mean <- (grandmean * foldchange) + rnorm(n = 1, mean = 0, sd = stddev.of.within.means)
            }
            temp.mean <- ifelse(temp.mean < 0, 0.0000001, temp.mean)
            temp.size <- exp(params$intraSD_dispersion.beta0 + (params$intraSD_dispersion.beta1 / temp.mean))
            if (distribution == "NB") {
                # sampled from NB
                temp.cells <- rnbinom(
                    n = design[which(design$cluster == cluster_id), "ncells"],
                    mu = temp.mean,
                    size = temp.size
                )
            } else if (distribution == "poisson") {
                # sampled from poisson
                temp.cells <- rpois(
                    n = design[which(design$cluster == cluster_id), "ncells"],
                    lambda = temp.mean
                )
            }

            temp.cells <- ifelse(
                rbinom(n = length(temp.cells), size = 1, prob = prob.zero) == 1,
                0,
                temp.cells
            )

            # FIXME: there is a dedicated function or at least a strategy for this but I forgot
            start <- 1 + sum(design$ncells[1:which(design$cluster == cluster_id) - 1])
            end <- sum(design$ncells[1:which(design$cluster == cluster_id)])

            simulatedData[i, start:end] <- temp.cells
        }
    }

    rownames(simulatedData) <- paste0("Gene_", seq_len(nrow(simulatedData)))
    colnames(simulatedData) <- paste0("Cell_", seq_len(ncol(simulatedData)))
    simulatedData <- as.matrix(simulatedData)

    simulatedData
}
