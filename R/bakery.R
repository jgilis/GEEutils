# Calculate the modified GEE variance estimator proposed by Kauermann and Carroll (2001).
#' @import geesmv
.bakery_KC <- function(gee.fit, data, mat, corstr) {

    # TODO: if this part turns out to be identical for all extraSandwich, consider moving it to a function .bakery_preprocess
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
        print(corstr)
        stop("'working correlation structure' not recognized")
    }
    #####

    cov.beta <- matrix(0, nrow = len, ncol = len)
    step11 <- matrix(0, nrow = len, ncol = len)
    for (i in 1:size) {
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(
            mat[, -length(mat[1, ])],
            mat$subj == unique(data$id)[i]
        ))
        var_i <- var[1:cluster$n[i], 1:cluster$n[i]]
        D <- mat.prod(covariate, exp(covariate %*% beta_est))
        Vi <- diag(
            sqrt(c(exp(covariate %*% beta_est))),
            cluster$n[i]
        ) %*% var_i %*% diag(sqrt(c(exp(covariate %*%
            beta_est))), cluster$n[i])
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
    }
    step12 <- matrix(0, nrow = len, ncol = len)
    step13 <- matrix(0, nrow = len_vec, ncol = 1)
    step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
    p <- matrix(0, nrow = len_vec, ncol = size)
    for (i in 1:size) {
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(
            mat[, -length(mat[1, ])],
            mat$subj == unique(data$id)[i]
        ))
        var_i <- var[1:cluster$n[i], 1:cluster$n[i]]
        D <- mat.prod(covariate, exp(covariate %*% beta_est))
        Vi <- diag(
            sqrt(c(exp(covariate %*% beta_est))),
            cluster$n[i]
        ) %*% var_i %*% diag(sqrt(c(exp(covariate %*%
            beta_est))), cluster$n[i])
        ##### KC #####
        xy <- t(D) %*% solve(Vi) %*% mat.sqrt.inv(cormax.ind(cluster$n[i]) -
            D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*%
            (y - exp(covariate %*% beta_est))
        ##############
        
        # check and handle complex numbers from eigenvalues
        if(any(is.complex(xy))){
            message(paste0("Subject ", i, ": Complex numbers detected"))
            xy <- zapsmall(xy)
            if (all(Im(xy)==0)){
                xy <- as.numeric(xy)
                message(paste0("Subject ", i, ": All complex +0i, set to numeric"))
            } else {
                message(paste0("Subject ", i, ": At least some complex not +0i"))
            }
        }
        
        step12 <- step12 + xy %*% t(xy)
        # TODO: could change matrixcalc::vec to c
        step13 <- step13 + matrixcalc::vec(xy %*% t(xy))
        p[, i] <- matrixcalc::vec(xy %*% t(xy))
    }
    for (i in 1:size) {
        dif <- (p[, i] - step13 / size) %*% t(p[, i] - step13 / size)
        step14 <- step14 + dif
    }
    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size / (size - 1) * kronecker(solve(step11), solve(step11)) %*%
        step14 %*% kronecker(solve(step11), solve(step11))

    return(cov.beta) # when not returning cov.var, some steps could be omitted
}


# Calculate the modified GEE variance estimator proposed by Pan (2001).
#' @import geesmv
.bakery_Pan <- function(gee.fit, data, mat, corstr) {

    # TODO: if this part turns out to be identical for all extraSandwich, consider moving it to a function .bakery_preprocess
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
        print(corstr)
        stop("'working correlation structure' not recognized")
    }
    #####

    cov.beta <- matrix(0, nrow = len, ncol = len)
    step <- matrix(0, nrow = cluster$n[1], ncol = cluster$n[1])
    for (i in 1:size) {
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(
            mat[, -length(mat[1, ])],
            mat$subj == unique(data$id)[i]
        ))
        resid <- (y - exp(covariate %*% beta_est)) %*% t(y -
            exp(covariate %*% beta_est))
        B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(B) <- 1 / sqrt(exp(covariate %*% beta_est))
        step <- step + B %*% resid %*% B
    }
    unstr <- step / size
    step11 <- matrix(0, nrow = len, ncol = len)
    step12 <- matrix(0, nrow = len, ncol = len)
    step13 <- matrix(0, nrow = len_vec, ncol = 1)
    step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
    p <- matrix(0, nrow = len_vec, ncol = size)

    for (i in 1:size) {
        y <- as.matrix(data$response[data$id == unique(data$id)[i]])
        covariate <- as.matrix(subset(
            mat[, -length(mat[1, ])],
            mat$subj == unique(data$id)[i]
        ))
        var_i <- var[1:cluster$n[i], 1:cluster$n[i]]

        A <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(A) <- exp(covariate %*% beta_est)
        D <- mat.prod(covariate, exp(covariate %*% beta_est))
        Vi <- diag(
            sqrt(c(exp(covariate %*% beta_est))),
            cluster$n[i]
        ) %*% var_i %*% diag(sqrt(c(exp(covariate %*%
            beta_est))), cluster$n[i])
        xy <- t(D) %*% solve(Vi) %*% sqrt(A) %*% unstr %*%
            sqrt(A) %*% solve(Vi) %*% D
        xx <- t(D) %*% solve(Vi) %*% D
        step12 <- step12 + xy
        step11 <- step11 + xx
        step13 <- step13 + matrixcalc::vec(xy)
        p[, i] <- matrixcalc::vec(xy)
    }

    for (i in 1:size) {
        dif <- (p[, i] - step13 / size) %*% t(p[, i] - step13 / size)
        step14 <- step14 + dif
    }

    cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
    cov.var <- size / (size - 1) * kronecker(solve(step11), solve(step11)) %*%
        step14 %*% kronecker(solve(step11), solve(step11))

    return(cov.beta) # when not returning cov.var, some steps could be omitted
}


#' Compute robust standard error estimates using multiple sandwich strategies
#'
#' @description Compute robust standard error estimates using multiple sandwich
#'   strategies
#'
#' @param formula A formula expression as for other regression models, of the
#'   form response ~ predictors. See the documentation of `lm` and `formula` for
#'   details.
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
#'   and variance functions. Currently, only "poisson" is supported.
#'
#' @param corstr A `character string` specifying the correlation structure.
#'   Default is "exchangeable". The following are permitted: "independence",
#'   "fixed", "stat_M_dep", "non_stat_M_dep", "exchangeable", "AR-M" and
#'   "unstructured".
#'
#' @param extraSandwich A `character vector` of sandwich estimation procedures
#'   that must be used to compute robust standard errors, additional to the
#'   model standard error and the robust estimator from Liang and Zeger. The
#'   following are permitted: "KC" and "Pan". If not specified or set to "none",
#'   only the model standard error and the robust estimator from Liang and Zeger
#'   will be provided.
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
#'     formula = gene ~ group_id,
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
#' @importFrom matrixcalc vec
#'
#' @export
bakery <- function(formula,
                   id, data,
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
        data$response <- model.response(m, "numeric")
        mat <- as.data.frame(model.matrix(formula, m))
    } else {
        m <- model.frame(formula, data)
        mt <- attr(m, "terms")
        data$response <- model.response(m, "numeric")
        mat <- as.data.frame(model.matrix(formula, m))
    }
    gee.fit <- gee(formula,
        data = data,
        id = id,
        family = family,
        corstr = corstr,
        silent = TRUE
    )

    if ("none" %in% extraSandwich) {
        return(gee.fit)
    }

    if ("KC" %in% extraSandwich) {
        gee.fit$KC.variance <- .bakery_KC(gee.fit, data, mat, corstr)
    }

    if ("Pan" %in% extraSandwich) {
        gee.fit$Pan.variance <- .bakery_Pan(gee.fit, data, mat, corstr)
    }
    # only return what will be used later
    gee.fit <- gee.fit[grepl("variance|coefficients", names(gee.fit))]
    return(gee.fit)
}
