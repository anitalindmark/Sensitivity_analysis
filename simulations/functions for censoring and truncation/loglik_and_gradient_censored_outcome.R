#' Implementation of loglikelihood function for ML estimation of regression parameters, censored outcome
#'
#'Implementation of loglikelihood function for ML estimation of regression parameters for a continuous mediator and a continuous, censored outcome. 
#'
#'@param par Vector of parameter values.
#'@param Rho The value of the sensitivity parameter.
#'@param X.expl The model matrix (see \code{\link{model.matrix}}) of the mediator model 
#'@param X.resp The model matrix (see \code{\link{model.matrix}}) of the outcome model
#'@param outc.resp The outcome, a vector.
#'@param outc.expl The mediator, a column matrix.
#'@param L Censoring point. Only left censoring is implemented in the function.
#'

LogL.cens <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl=outc.expl, L= L){
  
  # Separating the coefficients from the outcome (model.resp) and mediator (model.expl) models
  d.resp <- dim(X.resp)[2] 
  coef.resp <- as.matrix(par[1:d.resp])          # Outcome model parameters
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par) - 1)]) # Mediator model parameters
  sigma.res.resp <- par[d.resp + 1] # Error standard deviation outcome model
  sigma.res.expl <- par[length(par)] # Error standard deviation mediator model
  
  # Components to be used in the log-likelihood
  Q1 <- outc.expl - X.expl%*%expl.coef
  Q2 <- outc.resp - X.resp%*%coef.resp
  Qs <- cbind(Q1, Q2)
  
  Sigma <- cbind(c(sigma.res.expl^2, Rho*sigma.res.expl*sigma.res.resp), c(Rho*sigma.res.expl*sigma.res.resp, sigma.res.resp^2))
  
  # Calculating the log-likelihood:
  logl <- sum(ifelse(outc.resp > L, log(dmvnorm(Qs, sigma = Sigma)), log(pnorm((L - X.resp%*%coef.resp - 
               Rho*sigma.res.resp/sigma.res.expl*(outc.expl - X.expl%*%expl.coef))/(sigma.res.resp*sqrt(1 - Rho^2)))) - 
                 log(sigma.res.expl) + log(dnorm((outc.expl - X.expl%*%expl.coef)/sigma.res.expl))))
  
  return(logl) 
}

#' Analytic gradients of the loglikelihood function for ML estimation of regression parameters, censored outcome
#'
#'Implementation of the analytic gradients of the loglikelihood function for ML estimation of regression parameters for 
#'a continuous mediator and a continuous, censored outcome. 
#' 
#'@param par Vector of parameter values.
#'@param Rho The value of the sensitivity parameter.
#'@param X.expl The model matrix (see \code{\link{model.matrix}}) of the mediator model 
#'@param X.resp The model matrix (see \code{\link{model.matrix}}) of the outcome model
#'@param outc.resp The outcome, a vector.
#'@param outc.expl The mediator, a column matrix.
#'@param L Censoring point. Only left censoring is implemented in the function.
#'

grr.cens<- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, L = L){ 
  
  # Separating the coefficients from the outcome (model.resp) and mediator (model.expl) models
  d.resp<-dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp]) # Outcome model parameters         
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par) - 1)]) # Mediator model parameters
  sigma.res.resp <- par[d.resp + 1] # Error standard deviation outcome model
  sigma.res.expl <- par[length(par)] # Error standard deviation mediator model
  
  # Components to be used in the calculation of the gradients
  A <- (outc.expl - X.expl%*%expl.coef)/sigma.res.expl
  B <- (outc.resp - X.resp%*%coef.resp)/sigma.res.resp
  C <- (L - X.resp%*%coef.resp- Rho*sigma.res.resp/sigma.res.expl*(outc.expl - 
         X.expl%*%expl.coef))/(sigma.res.resp*sqrt(1 - Rho^2))
  
  n <- length(outc.resp)
   
  # Calculation of the gradients
  gr.expl.coef <- colSums(ifelse(outc.resp > L, 1, 0)*as.vector(A - Rho*B)*X.expl/((1 - Rho^2)*sigma.res.expl)) +
                    colSums(ifelse(outc.resp > L, 0, 1)*(as.vector(dnorm(C)/pnorm(C))*Rho*X.expl/(sigma.res.expl*sqrt(1 - 
                     Rho^2)) + as.vector(A/sigma.res.expl)*X.expl)) 
  
  gr.coef.resp <- colSums(ifelse(outc.resp > L, 1, 0)*as.vector(B - Rho*A)*X.resp/((1 - Rho^2)*sigma.res.resp)) -
                   colSums(ifelse(outc.resp > L, 0, 1)*as.vector(dnorm(C)/pnorm(C))*X.resp/(sigma.res.resp*sqrt(1 - Rho^2)))
                                                                                                                                                                                                                  
  gr.sigma.b <- -sum(ifelse(outc.resp > L, 1, 0))/sigma.res.expl + 
                  1/(sigma.res.expl^2*(1 - Rho^2))*sum(ifelse(outc.resp > L, 1, 0)*((outc.expl - 
                   X.expl%*%expl.coef)*(A - Rho*B))) + 
                    sum(ifelse(outc.resp > L, 0, 1)*(dnorm(C)/pnorm(C)*Rho*A/(sigma.res.expl*sqrt(1 - Rho^2)) + 
                     A^2/(sigma.res.expl) - 1/sigma.res.expl))
  
  gr.sigma.t <- -sum(ifelse(outc.resp > L, 1, 0))/sigma.res.resp + 1/(sigma.res.resp^2*(1 - 
                   Rho^2))*sum(ifelse(outc.resp > L, 1, 0)*((outc.resp - X.resp%*%coef.resp)*(B - Rho*A))) -
                     sum(ifelse(outc.resp > L, 0, 1)*dnorm(C)/pnorm(C)*(L - X.resp%*%coef.resp)/(sigma.res.resp^2*sqrt(1 - Rho^2)))

  return(c(gr.coef.resp, gr.sigma.t, gr.expl.coef, gr.sigma.b)) 
}

