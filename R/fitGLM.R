#' Fitting a Generalized Linear Model
#'
#' Fit a (quasi-)Poisson Generalized Linear Model using [stats::glm()], for
#' further use in DE testing.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object.
#' @param formula An object of class ["formula"][stats::formula] or a character
#'   string that can be coerced to one. Entries in the formula should refer to
#'   columns in `colData(sce)`.
#' @param family Character string specifying the distribution to be used for the
#'   GLM. Currently supports `"poisson"`, `"quasipoisson"` and `"negbin"`
#'   (negative binomial GLM, through [MASS::glm.nb()]).
#' @param offsets Either a character string specifying which method to use to
#'   calculate offsets or a numeric vector providing the offsets.
#' @param prior_count Prior count to be added to counts. Default: 0.
#' @param ... Further arguments passed on to [stats::glm()] or [MASS::glm.nb()]
#'   when `family = "negbin"`.
#'
#' @details
#' The `offsets` are calculated using [edgeR::calcNormFactors()] by default
#' (`method = "TMM"`). Alternatively, a vector of offsets can be provided with
#' length equal to the number of cells in the data.
#'
#' @return
#' A list of `"glm"` objects for each gene.
#'
#' @examples
#' ## Mock up data set
#' library(scuttle)
#' sce <- mockSCE(ncells = 100, ngenes = 500)
#' colData(sce)
#'
#' ## Fit model using available colData columns
#' model_fits <- fitGLM(sce, ~ Mutation_Status + Treatment)
#' class(model_fits)
#' length(model_fits)
#'
#' @export
#' @importFrom SingleCellExperiment counts
#' @importFrom stats as.formula
fitGLM <- function(sce, formula,
                   family = c("poisson", "quasipoisson", "negbin"),
                   offsets = "TMM", prior_count = 0,
                   ...) {

    stopifnot(is(sce, "SummarizedExperiment"))

    formula <- as.formula(formula)
    family <- match.arg(family)

    ## Requires raw counts
    Y <- counts(sce) + prior_count
    cd <- as.data.frame(colData(sce))

    offsets <- .handle_offsets(Y = Y, offsets = offsets)

    ## Fit GLM for each gene
    # TODO: this can be easily parallelized, look into beachmat::rowBlockApply()
    apply(Y, 1,
        FUN = .fit_glm,
        formula = formula,
        offsets = offsets,
        family = family,
        cd = cd,
        ...,
        simplify = FALSE
    )
}


## Helper to fit single GLM from count vector and formula
#' @importFrom stats glm
#' @importFrom MASS glm.nb
.fit_glm <- function(y, offsets, formula, family, cd, ...) {
    ## Make formula of the form "y ~ ..." and include offsets
    cd$offsets <- offsets
    cd$y <- y
    formula <- as.formula(
        paste("y", paste(formula, collapse = " "), "+ offset(offsets)")
    )

    if (identical(family, "negbin")) {
        ## Suppress partial argument match warnings
        suppressWarnings({
            out <- glm.nb(formula = formula, data = cd, ...)
        })
    } else {
        out <- glm(formula = formula, family = family, data = cd, ...)
    }
    out
}


.handle_offsets <- function(Y, offsets) {
    n_cells <- ncol(Y)

    if (is.null(offsets)) {
        offsets <- rep(0, n_cells)
    } else if (is.character(offsets) && length(offsets) == 1) {
        ## Currently only supports edgeR's norm.factors
        offsets <- .get_offsets(Y, method = offsets)
    } else if (is.numeric(offsets)) {
        if (!(length(offsets) == ncol(Y))) {
            stop(
                "Wrong length for `offsets`. Should be of length `ncol(sce)`.",
                call. = FALSE
            )
        }
    } else {
        stop(
            "`offsets` should be either a numeric vector or a character string.",
            call. = FALSE
        )
    }
    offsets
}


#' @importFrom Matrix colSums
#' @importFrom edgeR calcNormFactors
.get_offsets <- function(Y, method) {
    lib_sizes <- colSums(Y)
    norm_factors <- calcNormFactors(Y, lib.size = lib_sizes, method = method)
    log(lib_sizes) + log(norm_factors)
}
