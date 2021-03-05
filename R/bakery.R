# Calculate the modified GEE variance estimator proposed by 
# Kauermann and Carroll (2001). Code adapted from the `GEE.var.kc` function 
# of the geesmv package
#' @import geesmv
.bakery_KC <- function(gee.fit, data, mat, offset, corstr, family) {

    # TODO: if this part turns out to be identical for all extraSandwich, 
    # consider moving it to a function .bakery_preprocess
    beta_est <- gee.fit$coefficients
    alpha <- gee.fit$working.correlation[1, 2]
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- gee.fit$id
    cluster <- cluster.size(data$id)
    ncluster <- max(cluster$n)
    size <- cluster$m
    mat$subj <- rep(unique(data$id), cluster$n)
    
    if (is.character(corstr)) {
        var <- switch(corstr,
            independence = cormax.ind(ncluster),
            exchangeable = cormax.exch(ncluster, alpha),
            `AR-M` = cormax.ar1(ncluster, alpha),
            unstructured = summary(gee.fit)$working.correlation,
        )
    } else {
        stop("Correlation structure `", corstr, "` not recognized.", call. = FALSE)
    }
    
    if (is.character(family)) {
        family <- switch(family, poisson = "poisson", binomial = "binomial")
    }
    else {
        if (is.function(family)) {
            family <- family()[[1]]
        }
        else {
            stop("Family `", family, "` not recognized.", call. = FALSE)
        }
    }
    ###############

    cov.beta <- matrix(0, nrow = len, ncol = len)
    step11 <- matrix(0, nrow = len, ncol = len)
    
    for (i in 1:size) { #loop 1
        
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(mat[, -length(mat[1, ])], 
                                      mat$subj == unique(data$id)[i]))
        var_i <- var[1:cluster$n[i], 1:cluster$n[i]]
        
        off <- offset[mat$subj == unique(data$id)[i]]
        eta <- covariate %*% beta_est + off # add offsets (0 if not provided)
        
        if (family == "poisson") {
            mu <- exp(eta) # poisson link
            
            D <- mat.prod(covariate, mu)
            Vi <- diag(sqrt(c(mu)), cluster$n[i]) %*% 
                var_i %*% 
                diag(sqrt(c(mu)), cluster$n[i])
            
            xx <- t(D) %*% solve(Vi) %*% D
            step11 <- step11 + xx
            
        } else if (family == "binomial") {
            
            D <- mat.prod(covariate, exp(eta)/((1 + exp(eta))^2)) # binom link
            
            Vi <- diag(sqrt(c(exp(eta)/(1 + exp(eta))^2)), cluster$n[i]) %*% 
                var_i %*% 
                diag(sqrt(c(exp(eta)/(1 + exp(eta))^2)), cluster$n[i])
            
            xx <- t(D) %*% solve(Vi) %*% D
            step11 <- step11 + xx
        }
    }
    
    step12 <- matrix(0, nrow = len, ncol = len)
    # step13 <- matrix(0, nrow = len_vec, ncol = 1)
    # step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
    # p <- matrix(0, nrow = len_vec, ncol = size)
    
    # TODO: Several computations of the loop below are the same as in loop above
    #-> we can port the objects from above, especially solve(Vi), to loop 2. 
    # I did this locally -> ±2x speed gain
    for (i in 1:size) { #loop 2
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(mat[, -length(mat[1, ])], 
                                      mat$subj == unique(data$id)[i]))
        var_i <- var[1:cluster$n[i], 1:cluster$n[i]]
        
        off <- offset[mat$subj == unique(data$id)[i]]
        eta <- covariate %*% beta_est + off # add offsets (0 if not provided)
        
        if (family == "poisson") {
            mu <- exp(eta)
            
            D <- mat.prod(covariate, mu)
            Vi <- diag(sqrt(c(mu)), cluster$n[i]) %*% 
                var_i %*% 
                diag(sqrt(c(mu)), cluster$n[i])
            
            ######### core of KC algorithm #########
            xy <- t(D) %*% 
                solve(Vi) %*% 
                mat.sqrt.inv(cormax.ind(cluster$n[i]) - D %*% solve(step11) %*% 
                                 t(D) %*% solve(Vi)) %*% 
                (y - mu)
            ########################################
            
            # this part was added; not in geesmv
            if (any(is.complex(xy))) {
                message(paste0("Subject ", i, ": Complex numbers detected"))
                xy <- zapsmall(xy)
                if (all(Im(xy) == 0)) {
                    xy <- as.numeric(xy)
                    message(paste0("Subject ", i, ": All complex +0i, set to numeric"))
                }
                else {
                    message(paste0("Subject ", i, ": At least some complex not +0i"))
                }
            }
            step12 <- step12 + xy %*% t(xy)
            # step13 <- step13 + matrixcalc::vec(xy %*% t(xy))
            # p[, i] <- matrixcalc::vec(xy %*% t(xy))
            
        } else if (family == "binomial") {
            
            D <- mat.prod(covariate, exp(eta)/((1 + exp(eta))^2))
            Vi <- diag(sqrt(c(exp(eta)/(1 + exp(eta))^2)), cluster$n[i]) %*% 
                var_i %*% 
                diag(sqrt(c(exp(eta)/(1 + exp(eta))^2)), cluster$n[i])
            
            ######### core of KC algorithm #########
            xy <- t(D) %*% 
                solve(Vi) %*% 
                mat.sqrt.inv(cormax.ind(cluster$n[i]) - D %*% solve(step11) %*% 
                                 t(D) %*% solve(Vi)) %*% 
                (y - exp(eta)/(1 + exp(eta)))
            ########################################
            
            # this part was added; not in geesmv
            if (any(is.complex(xy))) {
                message(paste0("Subject ", i, ": Complex numbers detected"))
                xy <- zapsmall(xy)
                if (all(Im(xy) == 0)) {
                    xy <- as.numeric(xy)
                    message(paste0("Subject ", i, ": All complex +0i, set to numeric"))
                }
                else {
                    message(paste0("Subject ", i, ": At least some complex not +0i"))
                }
            }
            
            step12 <- step12 + xy %*% t(xy)
            # step13 <- step13 + matrixcalc::vec(xy %*% t(xy))
            # p[, i] <- matrixcalc::vec(xy %*% t(xy))
        }
    }
    
    # for (i in 1:size) {
    #     dif <- (p[, i] - step13 / size) %*% t(p[, i] - step13 / size)
    #     step14 <- step14 + dif
    # }
    
    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    # cov.var <- size / (size - 1) * kronecker(solve(step11), solve(step11)) %*%
    #     step14 %*% kronecker(solve(step11), solve(step11))

    return(cov.beta) # when not returning cov.var, some steps could be omitted
}


# Calculate the modified GEE variance estimator proposed by Pan (2001).
# Code adapted from the `GEE.var.pan` function of the geesmv package.
# Note that, in this implementation, all clusters must be the same size.
#' @import geesmv
.bakery_Pan <- function(gee.fit, data, mat, offset, corstr, family) {

    # TODO: if this part turns out to be identical for all extraSandwich, 
    # consider moving it to a function .bakery_preprocess
    beta_est <- gee.fit$coefficients
    alpha <- gee.fit$working.correlation[1, 2]
    len <- length(beta_est)
    len_vec <- len^2
    data$id <- gee.fit$id
    cluster <- cluster.size(data$id)
    ncluster <- max(cluster$n)
    size <- cluster$m
    mat$subj <- rep(unique(data$id), cluster$n)
    
    if (is.character(corstr)) {
        var <- switch(corstr,
            independence = cormax.ind(ncluster),
            exchangeable = cormax.exch(ncluster, alpha),
            `AR-M` = cormax.ar1(ncluster, alpha),
            unstructured = summary(gee.fit)$working.correlation,
        )
    } else {
        stop("Correlation structure `", corstr, "` not recognized.", call. = FALSE)
    }
    
    if (is.character(family)) {
        family <- switch(family, poisson = "poisson", binomial = "binomial")
    }
    else {
        if (is.function(family)) {
            family <- family()[[1]]
        }
        else {
            stop("Family `", family, "` not recognized.", call. = FALSE)
        }
    }
    ###############

    cov.beta <- matrix(0, nrow = len, ncol = len)
    step <- matrix(0, nrow = cluster$n[1], ncol = cluster$n[1])
    
    for (i in 1:size) { # loop 1
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(
            mat[, -length(mat[1, ])], 
            mat$subj == unique(data$id)[i]))
        
        off <- offset[mat$subj == unique(data$id)[i]]
        eta <- covariate %*% beta_est + off # add offsets (0 if not provided)
        
        if (family == "poisson") {
            mu <- exp(eta) # poisson link
            resid <- (y - mu) %*% t(y - mu)
            B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
            diag(B) <- 1/sqrt(mu)
            step <- step + B %*% resid %*% B
            
        } else if (family == "binomial") {
            resid <- (y - exp(eta)/(1 + exp(eta))) %*% 
                t(y - exp(eta)/(1 + exp(eta))) # binom link
            
            B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
            diag(B) <- 1/sqrt(exp(eta)/(1 + exp(eta))^2)
            step <- step + B %*% resid %*% B
        }
    }

    unstr <- step / size
    
    step11 <- matrix(0, nrow = len, ncol = len)
    step12 <- matrix(0, nrow = len, ncol = len)
    # step13 <- matrix(0, nrow = len_vec, ncol = 1)
    # step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
    # p <- matrix(0, nrow = len_vec, ncol = size)

    for (i in 1:size) {
        
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(mat[, -length(mat[1, ])], 
                                      mat$subj == unique(data$id)[i]))
        var_i <- var[1:cluster$n[i], 1:cluster$n[i]]
        
        off <- offset[mat$subj == unique(data$id)[i]]
        eta <- covariate %*% beta_est + off # add offsets (0 if not provided)
        
        if (family == "poisson") {
            mu <- exp(eta) # poisson link
            A <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
            diag(A) <- mu
            
            D <- mat.prod(covariate, mu)
            Vi <- diag(sqrt(c(mu)), cluster$n[i]) %*% 
                var_i %*% 
                diag(sqrt(c(mu)), cluster$n[i])
            
            xy <- t(D) %*% solve(Vi) %*% sqrt(A) %*% unstr %*% 
                sqrt(A) %*% solve(Vi) %*% D
            
            xx <- t(D) %*% solve(Vi) %*% D
            step12 <- step12 + xy
            step11 <- step11 + xx
            # step13 <- step13 + matrixcalc::vec(xy)
            # p[, i] <- matrixcalc::vec(xy)
            
        } else if (family == "binomial") {
            
            A <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
            diag(A) <- exp(eta)/(1 + exp(eta))^2 # binom link
            
            D <- mat.prod(covariate, exp(eta)/((1 + exp(eta))^2))
            Vi <- diag(sqrt(c(exp(eta)/(1 + exp(eta))^2)), cluster$n[i]) %*% 
                var_i %*% 
                diag(sqrt(c(exp(eta)/(1 + exp(eta))^2)), cluster$n[i])
            
            xy <- t(D) %*% solve(Vi) %*% sqrt(A) %*% unstr %*% 
                sqrt(A) %*% solve(Vi) %*% D
            
            xx <- t(D) %*% solve(Vi) %*% D
            step12 <- step12 + xy
            step11 <- step11 + xx
            # step13 <- step13 + matrixcalc::vec(xy)
            # p[, i] <- matrixcalc::vec(xy)
        }
    }

    # for (i in 1:size) {
    #     dif <- (p[, i] - step13 / size) %*% t(p[, i] - step13 / size)
    #     step14 <- step14 + dif
    # }

    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    # cov.var <- size / (size - 1) * kronecker(solve(step11), solve(step11)) %*%
    #     step14 %*% kronecker(solve(step11), solve(step11))

    return(cov.beta) # when not returning cov.var, some steps could be omitted
}


#' Compute robust standard error estimates using multiple sandwich strategies
#'
#' @description Compute robust standard error estimates using multiple sandwich
#'   strategies
#'   
#' @param response A `numeric` response vector. In the context of single-cell
#' transcriptomics, it is the expression vector of a single gene.
#'
#' @param formula A formula expression as for other regression models, of the
#'   form response ~ predictors. See the documentation of `lm` and `formula` for
#'   details. It is possible to provide offsets to the model by including 
#'   offset(my_offsets) to the model, with my_offsets a vector of precomputed
#'   offsets. 
#'
#' @param id A `vector` which identifies the clusters. The length of id should
#'   be the same as the number of observations. Data are assumed to be sorted so
#'   that observations on a cluster are contiguous rows for all entities in the
#'   formula.
#'
#' @param data An optional `data frame` in which to interpret the variables
#'   occurring in the formula, along with the id and n variables.
#'
#' @param family A `character string` indicating the family for defining link
#'   and variance functions. Currently, only "poisson" and "binomial are 
#'   supported.
#'
#' @param corstr A `character string` specifying the correlation structure.
#'   Default is "exchangeable". The following are permitted: "independence",
#'   "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and
#'   "unstructured".
#'
#' @param extraSandwich A `character vector` of sandwich estimation procedures
#'   that must be used to compute robust standard errors, additional to the
#'   model standard error and the robust estimator from Liang and Zeger. The
#'   following are permitted: "none", "DF", "KC" and "Pan". If not specified or 
#'   set to "none", only the naive model standard error and the robust estimator 
#'   from Liang and Zeger will be provided.
#'
#' @return An object of class "gee" including multiple sandwich estimators.
#'
#' @rdname bakery
#'
#' @author Jeroen Gilis
#'
#' @examples
#' ## Simulate single gene across 10 patients in 2 groups
#' n_cells <- 500
#' gene <- rpois(n_cells, 3)
#' group_id <- factor(rep(c("control", "treat"), each = n_cells / 2))
#' patient_id <- factor(rep(paste0("patient", 1:10), each = n_cells / 10))
#' data <- data.frame(gene, group_id, patient_id)
#'
#' ## Run bakery with KC and Pan extra sandwiches
#' out <- bakery(
#'     response = gene,
#'     formula = ~ group_id,
#'     id = "patient_id",
#'     data = data,
#'     family = "poisson",
#'     corstr = "exchangeable",
#'     extraSandwich = c("KC", "Pan")
#' )
#'
#' ## Liang-Zeger variance estimates
#' diag(out$robust.variance)
#'
#' ## KC variance estimates
#' diag(out$KC.variance)
#'
#' ## Pan variance estimates
#' diag(out$Pan.variance)
#'
#' @import gee
#'
#' @export
bakery <- function(response,
                   formula,
                   id, 
                   data,
                   family,
                   corstr,
                   extraSandwich = "none") {
    if (is.null(data$id)) {
        index <- which(names(data) == id)
        data$id <- data[, index]
    }
    init <- model.frame(formula, data)
    init$num <- 1:length(init[, 1])
    if (any(is.na(init))) {
        index <- na.omit(init)$num
        data <- data[index, ]
        m <- model.frame(formula, data)
        mt <- attr(m, "terms")
        data$response <- response
        mat <- as.data.frame(model.matrix(formula, m))
    } else {
        m <- model.frame(formula, data)
        mt <- attr(m, "terms")
        data$response <- response
        mat <- as.data.frame(model.matrix(formula, m))
    }
    
    nn <- dim(data)[1]
    off <- model.offset(m)
    if (is.null(off)) {
        offset <- rep(0, nn) 
        #if NULL, do not use offsets (set all to zero)
        #TODO alternatively, we could enforce a certain default offset
    } else {
        offset <- off
    }
    
    design <- model.matrix(formula, data = data)
    
    gee.fit <- gee(formula = response ~ -1 + design + offset(offset),
        data = data,
        id = id,
        family = family,
        corstr = corstr,
        silent = TRUE
    )

    if ("none" %in% extraSandwich) {
        return(gee.fit)
    }
    
    if ("DF" %in% extraSandwich) {
        K <- nlevels(as.factor(gee.fit$id)) # number of clusters
        p <- length(gee.fit$coefficients) # number of parameters
        gee.fit$DF.variance <- gee.fit$robust.variance * (K/(K-p))
    }

    if ("KC" %in% extraSandwich) {
        gee.fit$KC.variance <- .bakery_KC(gee.fit,
                                        data,
                                        mat,
                                        offset,
                                        corstr,
                                        family)
    }

    if ("Pan" %in% extraSandwich) {
        gee.fit$Pan.variance <- .bakery_Pan(gee.fit,
                                            data,
                                            mat,
                                            offset,
                                            corstr,
                                            family)
    }
    # only return what will be used later
    gee.fit <- gee.fit[grepl("variance|coefficients", names(gee.fit))]
    return(gee.fit)
}
