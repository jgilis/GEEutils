#' Helper function to sample (empirically) parameters from the observed
#' distribution of grandmeans and dropouts
#' TODO: theoretically, I would sample amount=nGenes. However, in some but not 
#' all cases, the approx function produces NA's, causing problems downstream.
#' A quick and dirty solution is to sample a much larger number (10K),
#' remove NAs, and sample nGenes from the remaining values.
#' In the future, we should check the origin of the NAs, but also discuss 
#' with Lieven if we should continue with the approx function (as suggested
#' by Lieven) or use a different strategy, e.g. bootstraps.
.sample_empirical <- function(sumStats, amount=10000){
  
  GrandMean <- sumStats[which(sumStats$GrandMean > 0),"GrandMean"]
  pdf_of_GrandMean <- density(GrandMean)
  
  sample_GrandMean <- approx(
    cumsum(pdf_of_GrandMean$y)/sum(pdf_of_GrandMean$y),
    pdf_of_GrandMean$x,
    runif(amount) # took the same number of values as in input data
  )$y
  
  sample_GrandMean <- sample_GrandMean[!is.na(sample_GrandMean)]
  sample_GrandMean <- sample_GrandMean[sample_GrandMean>1]
  sample_GrandMean <- sample(sample_GrandMean, nGenes)
  
  DropOut <- 1 - sumStats$DropOut
  DropOut <- DropOut[which(DropOut > 0)] # check from Zimmerman, necessary?
  
  pdf_of_DropOut <- density(DropOut)
  
  sample_Dropout <- approx(
    cumsum(pdf_of_DropOut$y)/sum(pdf_of_DropOut$y),
    pdf_of_DropOut$x,
    runif(nGenes) # took the same number of values as in input data
  )$y
  
  empirical_sample <- list(sample_GrandMean, sample_Dropout)
  names(empirical_sample) <- c("grandMean", "Dropout")
  
  return(empirical_sample)
}

#' Simulate a synthetic scRNA-Seq dataset based on parameters obtained with
#' getSimulationParameters a
#'
#' @description Simulate a synthetic scRNA-Seq dataset based on parameters
#'     obtained with getSimulationParameters and empirically obtained parameters
#'   
#' @param object A `SingleCellExperiment` instance.
#' 
#' @param clustering A `character vector` indicating the categorical variable
#'     that is used to color the dimensionality reduction plot. Must be one of
#'     `colnames(colData(object))`.
#'
#' @param design A `data frame` containing 3 columns. Column `cluster` is a
#'     `factor` that specifies the clustering of the observations. 
#'      Column `group` is a `factor` that specifies to which (treatment) group 
#'      the observations belong. Column `ncells` is a `numeric`indicating the 
#'      number of cells belonging to each of the clusters.
#'
#' @param SimulationParameters A `list` of simulation parameters obtained with
#'      the `GEEutils::getSimulationParameters`.
#'
#' @param nGenes A `numeric` indicating hoe many genes should be simulated.
#'
#' @param distribution A `character` indicating from which baseline distribution
#'      the data should be sampled. Either "NB" or "poisson".
#'
#' @param fractionDGE A `numeric` between 0 and 1, indicating the fraction of
#'      differentially expressed genes in the simulated data. Defaults to zero.
#' @param foldchange A `numeric` indicating the size of differential expression
#'      in terms of a fold change. Defaults to 1 (not differentially expressed,
#'      mock dataset).
#'
#' @return A `matrix` of simulated gene expression data.
#'
#' @rdname simulator
#'
#' @author Jeroen Gilis
#'
#' @examples
#' ## Original data to simulate from
#' sce <- scuttle::mockSCE(ncells = 100, ngenes = 1000)
#'
#' ## Add 4 random individuals
#' sce$patient_id <- factor(rep(paste0("patient", 1:4), each = ncol(sce) / 4))
#'
#' sim_params <- getSimulationParameters(
#'     sce, clustering = "patient_id",
#'     plot = FALSE, verbose = FALSE
#' )
#'
#' ## Set up design for simulator; 6 individuals across 2 groups
#' cluster <- factor(paste0("individual", 1:6))
#' group <- rep(c("group1", "group2"), each = length(cluster) / 2)
#'
#' ## Set total number of desired cells for each individual
#' ncells <- rep(50, length(cluster))
#'
#' design <- data.frame(cluster, group, ncells)
#'
#' ## Simulate using Poisson model
#' sim_pois <- empirical_simulator(
#'     object = sce,
#'     clustering = "patient_id",
#'     design = design,
#'     SimulationParameters = sim_params,
#'     nGenes = nrow(sce),
#'     distribution = "poisson",
#'     fractionDGE = 0.1, #' fraction DE genes between groups
#'     foldchange = 1.5
#' )
#'
#'
#'
#' @import stats
#'
#' @export
empirical_simulator <- function(object,
                                clustering,
                                design, 
                                SimulationParameters,
                                nGenes, 
                                distribution,
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
  simulatedData <- as.data.frame(matrix(data = NA, nrow = nGenes, 
                                        sum(design$ncells)))
  
  sumStats <- .getSumStats(object = object,
                           clustering = clustering)
  
  ## sample empirically
  empirical_sample <- .sample_empirical(sumStats)
  
  for (i in seq_len(nGenes)) {
    
    # empirically get grandmean and relationship with inter-individual standard deviations
    grandmean <- empirical_sample$grandMean[i]
    stddev.of.within.means <- SimulationParameters$grandmean_interSD.beta1 * 
      grandmean
    
    if (i <= nGenes * fractionDGE) {# its a DE gene
      # Now all foldchanges will be "foldchange" or "-foldchange".
      # TODO: It would be better to sample foldchanges from a distribution with mean "foldchange"
      foldchange <- ifelse(stats::rbinom(n = 1, size = 1, 
                                         prob = 0.5) == 1, foldchange, 1/foldchange)
    }
    else {
      foldchange <- 1
    }
    
    prob.zero <- empirical_sample$Dropout[i]
    
    prob.zero <- ifelse(prob.zero > 1, 
                        rgamma(n = 1, 
                               shape = SimulationParameters$drop.shape, 
                               rate = SimulationParameters$drop.rate), 
                        prob.zero)
    
    drop.sd <- SimulationParameters$drop.beta0 + SimulationParameters$drop.beta1 * 
      prob.zero + SimulationParameters$drop.beta2 * (prob.zero^2)
    drop.sd <- ifelse(drop.sd < 0, 0, drop.sd)
    prob.zero <- rnorm(n = 1, mean = prob.zero, sd = drop.sd)
    prob.zero <- ifelse(prob.zero < 0, 0, prob.zero)
    prob.zero <- ifelse(prob.zero > 1, 1, prob.zero)
    prob.zero <- 1 - prob.zero
    
    for (cluster_id in levels(design$cluster)) {
      reference <- design$group[1]
      if (design[which(design$cluster == cluster_id), "group"] == reference) {
        temp.mean <- grandmean + rnorm(n = 1, mean = 0, 
                                       sd = stddev.of.within.means)
      }
      else {
        temp.mean <- (grandmean * foldchange) + rnorm(n = 1, 
                                                      mean = 0, sd = stddev.of.within.means)
      }
      temp.mean <- ifelse(temp.mean < 0, 1e-07, temp.mean)
      temp.size <- exp(SimulationParameters$intraSD_dispersion.beta0 + 
                         (SimulationParameters$intraSD_dispersion.beta1/temp.mean))
      
      if (distribution == "NB") {
        temp.cells <- rnbinom(n = design[which(design$cluster == 
                                                 cluster_id), "ncells"], mu = temp.mean, size = temp.size)
      }
      else if (distribution == "poisson") {
        temp.cells <- rpois(n = design[which(design$cluster == 
                                               cluster_id), "ncells"], lambda = temp.mean)
      }
      
      temp.cells <- ifelse(rbinom(n = length(temp.cells), 
                                  size = 1, 
                                  prob = prob.zero) == 1, 0, 
                           temp.cells)
      # FIXME: there is a dedicated function or at least a strategy for this but I forgot
      start <- 1 + sum(design$ncells[1:which(design$cluster == 
                                               cluster_id) - 1])
      end <- sum(design$ncells[1:which(design$cluster == 
                                         cluster_id)])
      simulatedData[i, start:end] <- temp.cells
    }
  }
  
  rownames(simulatedData) <- paste0("Gene_", 1:nrow(simulatedData))
  colnames(simulatedData) <- paste0("Cell_", 1:ncol(simulatedData))
  simulatedData <- as.matrix(simulatedData)
  
  simulatedData
}
