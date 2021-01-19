# Compute summary statistics of the scRNA-Seq dataset that will later be used to compute simulation parameters
.getSumStats <- function(object, clustering) {
    # FIXME: use of `@` heavily discouraged; use appropriate accessor instead
    counts <- as.data.frame(as.matrix(object@assays@data[[1]]))
    nclust <- length(levels(as.factor(object[[clustering]])))

    computevar <- function(x) {
        tapply(x, object[[clustering]], function(a) {
            stats::var(a[a != 0])
        })
    }
    temp.intravar <- apply(counts, 1, computevar)

    computemeans <- function(x) {
        tapply(x, object[[clustering]], function(a) {
            mean(a[a != 0])
        })
    }
    temp.intrameans <- apply(counts, 1, computemeans)

    computedrop <- function(x) {
        tapply(x, object[[clustering]], function(a) {
            length(a[which(a == 0)]) / length(a)
        })
    }
    temp.drop <- apply(counts, 1, computedrop)

    temp.disp <- (temp.intrameans**2) / (temp.intravar - temp.intrameans) # computes dispersion as mean squared divided by var - mean

    sumStats <- as.data.frame(
        t(rbind(
            temp.intrameans, temp.intravar, temp.drop
        ))
    )
    colnames(sumStats) <- paste0(
        rep(c("Mean_", "Var_", "Drop_"), each = nclust),
        colnames(sumStats)
    )

    sumStats$InterStD <- apply(sumStats[, 1:nclust], 1, function(a) {
        stats::sd(a[!is.na(a)])
    })
    sumStats$GrandMean <- apply(sumStats[, 1:nclust], 1, function(a) {
        mean(a[!is.na(a)])
    })
    sumStats$IntraVar <- apply(
      sumStats[, (1 + nclust):(nclust * 2)], 1, median
    )
    sumStats$DropOut <- apply(
      sumStats[, (1 + (nclust * 2)):(nclust * 3)], 1, mean
    )
    sumStats$DropOutStD <- apply(
      sumStats[, (1 + (nclust * 2)):(nclust * 3)], 1, function(a) {
        stats::sd(a[!is.na(a)])
      }
    )
    sumStats$Dispersion <- apply(temp.disp, 2, mean)

    sumStats
}

# Compute hyperparameters related to dropout
.getHyperDropout <- function(sumStats, plot) {

    # Dropout mean
    sumStats$DropOut <- 1 - sumStats$DropOut
    sumStats <- sumStats[which(sumStats$GrandMean > 0), ]
    sumStats <- sumStats[which(sumStats$DropOut > 0), ]

    temp.fit.gamma <- suppressWarnings(
        fitdistrplus::fitdist(sumStats$DropOut, "gamma", method = "mle")
    )
    drop.shape <- temp.fit.gamma$estimate[1]
    drop.rate <- temp.fit.gamma$estimate[2]

    if (plot) {
        gg <- ggplot(sumStats, aes(DropOut)) +
          geom_histogram(
            aes(y = after_stat(density)),
            fill = "cornflowerblue", color = "black"
          ) +
          ylab("Density") +
          stat_function(
            fun = dgamma, args = list(shape = drop.shape, rate = drop.rate)
          ) +
          theme_bw()
          # TODO: move plotting to separate function?
          # cfr. https://r-pkgs.org/r.html?q=side%20effects#isolate-side-effects
          print(gg)
    }

    # Dropout StD
    temp.mean.var <- glm(DropOutStD ~ DropOut + I(DropOut**2),
        family = gaussian(link = "identity"),
        data = sumStats
    )
    temp.model.summary <- summary(temp.mean.var)
    dropoutstd.beta0 <- temp.model.summary$coefficients[1, 1]
    dropoutstd.beta1 <- temp.model.summary$coefficients[2, 1]
    dropoutstd.beta2 <- temp.model.summary$coefficients[3, 1]

    dropoutParams <- list(
        unname(drop.shape), unname(drop.rate),
        dropoutstd.beta0, dropoutstd.beta1, dropoutstd.beta2
    )
    names(dropoutParams) <- c(
        "drop.shape", "drop.rate",
        "drop.beta0", "drop.beta1", "drop.beta2"
    )

    dropoutParams
}

# Compute hyperparameters related to mean expression levels over all cells
.getHyperGrandmean <- function(sumStats, plot) {
    sumStats <- sumStats[which(sumStats$GrandMean > 0), ]

    temp.fit.gamma <- suppressWarnings(
        fitdistrplus::fitdist(sumStats$GrandMean, "gamma", method = "mle")
    )
    grandmean.shape <- temp.fit.gamma$estimate[1]
    grandmean.rate <- temp.fit.gamma$estimate[2]

    if (plot) {
        gg <- ggplot(sumStats, aes(GrandMean)) +
          geom_histogram(
              aes(y = after_stat(density)),
              fill = "cornflowerblue", color = "black"
          ) +
          ylab("Density") +
          stat_function(
              fun = dgamma,
              args = list(shape = grandmean.shape, rate = grandmean.rate)
          ) +
          ggtitle("Grandmean") +
          theme_bw()
        print(gg)

        gg <- ggplot(sumStats, aes(GrandMean)) +
            geom_histogram(
                aes(y = after_stat(density)),
                fill = "cornflowerblue", color = "black"
            ) +
            stat_function(
                fun = dgamma,
                args = list(shape = grandmean.shape, rate = grandmean.rate)
            ) +
            xlim(0, quantile(sumStats$GrandMean, 0.95)) +
            ylab("Density") +
            ggtitle("Grandmean (max 95% quant)") +
            theme_bw()
        print(gg)
    }

    grandmeanParams <- list(unname(grandmean.shape), unname(grandmean.rate))
    names(grandmeanParams) <- c("grandmean.shape", "grandmean.rate")

    grandmeanParams
}

# Compute hyperparameters describing the relationship between mean expression and inter-individual variability
.getHyperRelationship1 <- function(sumStats, plot) {
    temp.mean.inter <- glm(
      data = sumStats, InterStD ~ 0 + GrandMean,
      family = gaussian(link = "identity")
    )
    inter.beta1 <- summary(temp.mean.inter)$coefficients[1, 1]

    if (plot) {
        gg <- ggplot(sumStats, aes(x = GrandMean, y = InterStD)) +
            geom_point() +
            stat_smooth(
                method = "glm", formula = y ~ 0 + x,
                method.args = list(family = gaussian(link = "identity")),
                size = 1.6
            ) +
          ggtitle("GrandMean vs. InterStD") +
          theme_bw()
        print(gg)

        gg <- ggplot(sumStats, aes(x = GrandMean, y = InterStD)) +
            geom_point() +
            stat_smooth(
                method = "glm", formula = y ~ 0 + x,
                method.args = list(family = gaussian(link = "identity")),
                size = 1.6
            ) +
            xlim(0, quantile(sumStats$GrandMean, 0.95)) +
            ylim(0, quantile(sumStats$InterStD, 0.95)) +
            ggtitle("GrandMean vs. InterStD (both max 95% quant)") +
            theme_bw()
        print(gg)
    }

    relationship1Params <- list(inter.beta1)
    names(relationship1Params) <- c("grandmean_interSD.beta1")

    relationship1Params
}


# Compute hyperparameters describing the relationship between intra-individual variability and expression dispersion
.getHyperRelationship2 <- function(sumStats, plot) {
    intravar <- na.omit(as.data.frame(cbind(
        unlist(sumStats[, which(grepl("Mean_", colnames(sumStats)))]),
        unlist(sumStats[, which(grepl("Var_", colnames(sumStats)))])
    )))
    colnames(intravar) <- c("IntraMean", "IntraVar")
    intravar$Dispersion <- (intravar$IntraMean**2) / ((intravar$IntraVar) - intravar$IntraMean)

    intravar <- intravar[which(intravar$Dispersion > 0 & intravar$Dispersion < 1000), ]
    # This is a very important step for the Kang dataset.
    # Only very few observations (±15%) have an intravar that satisfies these criteria.
    # This is probably due to inherent diferences between smartseq2 and 10X data.

    temp.mean.intra <- glm(
        data = intravar,
        Dispersion ~ I(1 / IntraMean),
        family = gaussian(link = "log")
    )
    intra.beta0 <- summary(temp.mean.intra)$coefficients[1, 1]
    intra.beta1 <- summary(temp.mean.intra)$coefficients[2, 1]

    if (plot) {
        gg <- ggplot(intravar, aes(x = IntraMean, y = Dispersion)) +
            geom_point() +
            stat_smooth(
                method = "glm", formula = y ~ I(1 / x),
                method.args = list(family = gaussian(link = "log")),
                size = 1.6
            ) +
            ggtitle("IntraMean vs. dispersion") +
            theme_bw()
        plot(gg)

        gg <- ggplot(intravar, aes(x = IntraMean, y = Dispersion)) +
            geom_point() +
            stat_smooth(
                method = "glm",
                formula = y ~ I(1 / x),
                method.args = list(family = gaussian(link = "log")),
                size = 1.6
            ) +
            xlim(1, quantile(intravar$IntraMean, 0.98)) +
            ylim(0, quantile(intravar$Dispersion, 0.98)) +
            ggtitle("IntraMean vs. dispersion (both 98% quant)") +
            theme_bw()
        plot(gg)
    }

    relationship2Params <- list(intra.beta0, intra.beta1)
    names(relationship2Params) <- c("intraSD_dispersion.beta0", "intraSD_dispersion.beta1")

    relationship2Params
}

#' Obtain parameters for simulating synthetic asRNA-Seq data based on a real
#' scRNA-Seq dataset.
#'
#' @description Extract parameters from a scRNA-Seq dataset that are relevant
#'   for simulating synthetic asRNA-Seq data.
#'
#' @param object A `SingleCellExperiment` instance.
#'
#' @param clustering A `character vector` indicating the categorical variable
#'   that is used to color the dimensionality reduction plot. Must be one of
#'   `colnames(colData(object))`.
#'
#' @param plot `Logical` TRUE or FALSE; should intermediate visualizations be
#'   generated?
#'
#' @param verbose `Logical` TRUE or FALSE; should progress be printed?
#'
#' @return A `list` of 10 parameters that can be used to simulate synthetic
#'   scRNA-Seq datasets
#'
#' @rdname getSimulationParameters
#'
#' @author Jeroen Gilis
#'
#' @import ggplot2
#' @importFrom fitdistrplus fitdist
#'
#' @export
getSimulationParameters <- function(object,
                                    clustering,
                                    plot = TRUE,
                                    verbose = FALSE) {
    if (verbose) message("Summary statistics")
    sumStats <- .getSumStats(object, clustering)

    if (verbose) message("Hyperparameters dropout")
    dropoutParams <- .getHyperDropout(sumStats, plot)

    if (verbose) message("Hyperparameters grand mean")
    grandmeanParams <- .getHyperGrandmean(sumStats, plot)

    if (verbose) message("Hyperparameters relationship 1")
    relationship1Params <- .getHyperRelationship1(sumStats, plot)

    if (verbose) message("Hyperparameters relationship 2")
    relationship2Params <- .getHyperRelationship2(sumStats, plot)

    if (verbose) message("Combine all hyperparameters")
    SimulationParameters <- c(
        dropoutParams, grandmeanParams,
        relationship1Params, relationship2Params
    )

    SimulationParameters
}
