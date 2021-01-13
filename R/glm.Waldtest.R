# Extract model coefficients from a glm
.getCoef <- function(model){
  return(model$coefficients)
}

# Compute estimates for target contrast
.getEstimates <- function(model,contrast){
  coef <- .getCoef(model)
  if (is.null(coef)){
    coef <- rep(NA, times=nrow(contrast))
  }
  return(contrast%*%coef)
}

# Compute the standard error on a target model estimate
.varContrast <- function(model,contrast){
  if (!any(class(model)=="fitError")){
    sx <- summary(model) # maybe I can also only compute vcov to gain speed
    if (nrow(sx$cov.scaled) == length(contrast)){
      return(sqrt(diag(t(contrast) %*% sx$cov.scaled %*% contrast))) # scaled by dispersion
    }
  }
  return(NA)
}

# Compute the robust standard error on a target model estimate
.sandwichVarContrast <- function(model,contrast,adjust){
  if (!any(class(model)=="fitError")){
    if(adjust){
      sanwich_var <- sandwich(model, adjust=TRUE) # small sample adjustment, divide by n-k
    } else {
      sanwich_var <- sandwich(model, adjust=FALSE)
    }
    if (nrow(sanwich_var) == length(contrast)){
      return(sqrt(diag(t(contrast) %*% sanwich_var %*% contrast)))
    }
  }
  return(NA)
}

# Extract the degrees of freedom from a glm
.getDf <- function(model){
  return(model$df.residual)
}


#' Perform a Wald test on a glm with the possible of adopting robust standard errors
#'
#' @description Perform a Wald test on a glm with the possible of adopting robust standard errors
#'
#' @param models A list of generalized linear model. If a single model, make sure to input it as a list by `list(glm)`.
#'
#' @param contrast A `matrix` specifying the contrast, i.e. combinations of model parameters, of interest.
#'
#' @param sandwich A `logical` variable; should robust standard errors be computed (sandwich procedure by Liang and Zeger).
#'  Default is TRUE.
#'
#' @param adjust A `logical` variable; should a finite sample adjustment be made?
#'  This amounts to multiplication with n/(n-k) where n is the number of observations and k the number of estimated parameters.
#'
#' @examples TO DO
#'
#' @return A `Dataframe` containing the requested model parameter estimates, (robust) standard errors, degrees of freedom,
#'  Wald test statistics and p-values.
#'
#' @rdname glm.Waldtest
#'
#' @author Jeroen Gilis
#'
#' @export
glm.Waldtest <- function(models, contrast, sandwich = TRUE, adjust = FALSE){

  estimates <- sapply(models, .getEstimates, contrast = contrast)
  if(sandwich){
      se <- sapply(models, .sandwichVarContrast ,contrast = contrast, adjust) # uses sandwich SEs
  } else {
      se <- sapply(models, .varContrast, contrast = contrast) # uses GLM SEs
  }
  dfs <- sapply(models, .getDf)

  W_stats <- estimates / sqrt(se) # Wald statistics
  pvals <- 2 * pt(abs(W_stats), df = dfs, lower.tail = FALSE)

  return(data.frame(estimates, se, dfs, W_stats, pvals))
}
