#'Function for ML estimation of regression parameters for sensitivity analysis, censored outcome
#'
#'Function for ML estimation of regression parameters for sensitivity analysis for a continuous mediator and a continuous, censored outcome. 
#'The optimization is performed using \code{\link{maxLik}} and the "BFGS" method. 
#'
#'@param model.expl Fitted \code{\link{glm}} mediator model object with \code{family = "gaussian"}, and \code{link = "identity"}.
#'@param model.resp Fitted outcome model object from the \code{tobit} function of the \code{AER] package with \code{dist = "gaussian"}. 
#'@param Rho The sensitivity parameter vector, the correlation between the error terms in the mediator and outcome models. 
#'@param L Censoring point. Only left censoring is implemented in the function.
#'@param progress Logical, indicating whether or not the progress (i.e. the \code{\link{proc.time}} for each \code{Rho}) of the optimization will be output
#'@param ... Additional arguments to be passed on to the \code{maxLik} function. Can be used to set the \code{method} and \code{control} arguments of the \code{maxLik} function.
#' 

ML.cens <- function(model.expl, model.resp, Rho, L = 0, progress = TRUE, ...){ 
  
  X.expl <- as.matrix(stats::model.matrix(model.expl))
  X.resp <- as.matrix(stats::model.matrix(model.resp)) 
  
  outc.resp <- model.resp$y[,1]
  outc.expl <- as.matrix(model.expl$y)
  
  
  nrho <- length(Rho)
  p1 <- length(coef(model.resp)) + 1
  p2 <- length(coef(model.expl)) + 1
  
  coefs <- matrix(nrow = (p1  + p2), ncol = nrho) # vector with all model parameters
  rownames(coefs) <- c(names(coef(model.resp)), "sigma.res.resp", names(model.expl$coef), "sigma.res.expl")
  
  dots <- list(...)
  if(is.null(dots$method)) # which optimization method should be used?
    method <- "BFGS" #"NR"
  else
    method <- dots$method
  if(is.null(dots$control)) # any control arguments given?
    control <- NULL
  else
    control <- dots$control
  
  i0 <- which(Rho == 0)
  coefs[, i0] <- c(coef(model.resp), model.resp$scale, coef(model.expl), 
                   sqrt(summary(model.expl)$dispersion)) # storing the parameters for Rho = 0

  # Convergence information, matrix and stored values for Rho = 0
  max.info <- list(converged = array(dim=nrho), message = array(dim=nrho), method = array(dim=nrho),
                   iterations = array(dim=nrho))
  max.info$converged[i0] <- " "
  max.info$message[i0] <- " "
  max.info$method[i0] <- " "
  max.info$iterations[i0] <- NA
  
  value <- vector(length = nrho) # vector to store values of the log-likelihood
  names(value) <- paste(Rho)
  value[i0] <- LogL.cens(coefs[, i0], Rho = 0, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp,
                       outc.expl = outc.expl, L = L) # Value of the log-likelihood for Rho = 0
  
  sigma <- list() # List to store the covariance matrices for the model parameters
  
  # Storing the covariance matrices for Rho = 0:
  sigma[[i0]] <- matrix(0, nrow = p1 + p2, ncol = p1 + p2)
  glmVcov <- stats::vcov(model.resp)
  expl.glmVcov <- stats::vcov(model.expl)
  sigma[[i0]][1:(p1 ), 1:(p1)] <- glmVcov
  sigma[[i0]][(p1 + 1):(p1 + p2 - 1), (p1 + 1):(p1 + p2 - 1)] <- expl.glmVcov
  sigma[[i0]][(p1 + p2), (p1 + p2)] <- summary(model.expl)$dispersion/(2*model.expl$df.residual)
  dimnames.sigma <- c(dimnames(glmVcov)[[1]], dimnames(expl.glmVcov)[[1]], "sigma.res.expl") 
  dimnames(sigma[[i0]]) <- list(dimnames.sigma, dimnames.sigma)

  # Optimization for Rho < 0:
  if(sum(Rho < 0) > 0){
    for(i in 1:(i0-1)){
      
      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i0 - i], sep=" ", "\n")
        ptm <- proc.time()
      }
      
      # Functions for use by maxLik, log-likelihood and analytic gradient
      f <- function(par)
        LogL.cens(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, L = L)
      
      g <- function(par)
        grr.cens(par, Rho = Rho[i0 - i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, L = L)
      
           
      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, start = coefs[, i0 - i + 1], method = method)
      
      # Storing information from maximization
      max.info$converged[i0 - i] <- ML$code <= 2
      max.info$message[i0 - i] <- ML$message
      max.info$method[i0 - i] <- ML$type
      max.info$iterations[i0 - i] <- ML$iterations
      
      # Storing parameters and log-likelihood value from maximization
      coefs[,i0 - i] <- ML$estimate
      value[i0 - i] <- ML$maximum
      
      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i0 - i]] <- -solve(ML$hessian)
      
      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }  
  
  
  # Optimization for Rho > 0:
  if(sum(Rho > 0) > 0){
    for(i in (i0 + 1):nrho){
      
      if(progress == TRUE){
        cat("Optimization for Rho =", Rho[i], sep=" ", "\n")
        ptm <- proc.time()
      }
      
      # Functions for use by maxLik, log-likelihood and analytic gradient
      f <- function(par)
        LogL.cens(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, L = L)
      
      g <- function(par)
        grr.cens(par, Rho = Rho[i], X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, L = L)
      
      
      # Maximization
      ML <- maxLik::maxLik(logLik = f, grad = g, start = coefs[, i - 1], method = method)
      
      # Storing information from maximization
      max.info$converged[i] <- ML$code <= 2
      max.info$message[i] <- ML$message
      max.info$method[i] <- ML$type
      max.info$iterations[i] <- ML$iterations
      
      # Storing parameters and log-likelihood value from maximization
      coefs[, i] <- ML$estimate
      value[i]  <-  ML$maximum
      
      # Solving the hessian to obtain covariance matrix for the parameters
      sigma[[i]] <- -solve(ML$hessian)
      
      if(progress == TRUE)
        cat("   Time elapsed:", (proc.time()-ptm)[3], "s", "\n")
    }
  }
  
  # Storing the estimated parameters for model.expl and model.resp over Rho.
  expl.coef <-as.matrix(coefs[(p1 + 1):(p1 + p2 - 1), ])
  if(p2 == 2){
    expl.coef  <- t(expl.coef)
    rownames(expl.coef) <- names(model.expl$coefficients)
  }
  
  sigma.res.expl <- as.matrix(coefs[p1 + p2, ])
  
  coef <- as.matrix(coefs[1:(p1 - 1), ])
  if(p1 == 2){
    coef <- t(coef)
    rownames(coef) <- names(coef(model.resp))
  }
  
  sigma.res.resp <- as.matrix(coefs[p1, ])
  
  colnames(expl.coef) <- rownames(sigma.res.expl) <- rownames(sigma.res.resp) <- colnames(coef) <- names(sigma) <- paste(Rho)
  
  max.info <- lapply(max.info, stats::setNames, nm = paste(Rho))
  
  # Setting family and link arguments of the outcome model to enable use of the calc.effects function from the
  # sensmediation package to estimate NIE, NDE and standard errors
  model.resp$family <- list()
  model.resp$family$family <- "gaussian"
  model.resp$family$link <- "identity"
  
  return(list(coef = coef, Rho = Rho, expl.coef = expl.coef, model.expl = model.expl,
              model.resp = model.resp, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp,
              outc.expl = outc.expl, value = value, sigma.res.resp = sigma.res.resp, sigma.res.expl = sigma.res.expl,
              sigmas = sigma, max.info = max.info))
  
}


