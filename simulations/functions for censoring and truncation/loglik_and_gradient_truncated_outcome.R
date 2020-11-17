#' Implementation of loglikelihood function for ML estimation of regression parameters, truncated outcome
#'
#'Implementation of loglikelihood function for ML estimation of regression parameters a continuous mediator and a continuous, truncated outcome. 
#'
#'@param par Vector of parameter values.
#'@param Rho The value of the sensitivity parameter.
#'@param X.expl The model matrix (see \code{\link{model.matrix}}) of \code{model.expl}
#'@param X.resp The model matrix (see \code{\link{model.matrix}}) of \code{model.resp}
#'@param outc.resp The outcome of \code{model.resp}, a vector.
#'@param outc.expl The outcome of \code{model.expl}, a column matrix.
#'@param L Truncation point. Only left censoring is implemented in the function.
#'@param med.name A character string indicating the name of the mediator used in the models.
#'

LogL.trunc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, L = L,
                       med.name = med.name){
  
  # Separating the coefficients from the outcome (model.resp) and mediator (model.expl) models
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp])   # Outcome model parameters       
  coef.resp.wo2<-as.matrix(coef.resp[-which(names(par)==med.name)]) # Outcome model parameters except delta2
  coef.resp.th2 <- par[med.name] # delta2
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par) - 1)]) # Mediator model parameters
  sigma.resp <- par[d.resp + 1] # Error standard deviation outcome model
  sigma.expl <- par[length(par)] # Error standard deviation mediator model
  
  # Components to be used in the log-likelihood  
  X.resp.wo2 <- X.resp[,-which(names(par)==med.name)]
  Q1 <- outc.expl-X.expl%*%expl.coef
  Q2 <- outc.resp-X.resp%*%coef.resp
  Qs <- cbind(Q1,Q2)
  
  Sigma <- cbind(c(sigma.expl^2,Rho*sigma.expl*sigma.resp), c(Rho*sigma.expl*sigma.resp,sigma.resp^2))
  
  # Calculating the log-likelihood:
  logl <- sum(log(dmvnorm(Qs, sigma = Sigma)) -
                log(1 - pnorm((L - (X.resp.wo2%*%coef.resp.wo2 +
                 coef.resp.th2*(X.expl%*%expl.coef)))/sqrt(coef.resp.th2^2*sigma.expl^2 +
                  sigma.resp^2 + 2*coef.resp.th2*Rho*sigma.expl*sigma.resp))))
  
  return(logl)
}

#' Analytic gradients of the loglikelihood function for ML estimation of regression parameters, truncated outcome
#'
#'Implementation of the analytic gradients of the loglikelihood function for ML estimation of regression parameters for 
#'a continuous mediator and a continuous, truncated outcome. 
#' 
#'@param par Vector of parameter values.
#'@param Rho The value of the sensitivity parameter.
#'@param X.expl The model matrix (see \code{\link{model.matrix}}) of the mediator model 
#'@param X.resp The model matrix (see \code{\link{model.matrix}}) of the outcome model
#'@param outc.resp The outcome, a vector.
#'@param outc.expl The mediator, a column matrix.
#'@param L Truncation point. Only left truncation is implemented in the function.
#'@param med.name A character string indicating the name of the mediator used in the models.
#'

grr.trunc <- function(par, Rho, X.expl = X.expl, X.resp = X.resp, outc.resp = outc.resp, outc.expl = outc.expl, L = L,
                      med.name = med.name){ 
  
  # Separating the coefficients from the outcome (model.resp) and mediator (model.expl) models
  d.resp <- dim(X.resp)[2]
  coef.resp <- as.matrix(par[1:d.resp]) # Outcome model parameters           
  coef.resp.wo2<-as.matrix(coef.resp[-which(names(par)==med.name)]) # Outcome model parameters except delta2
  coef.resp.th2 <- par[med.name] # delta2
  expl.coef <- as.matrix(par[(d.resp + 2):(length(par) - 1)]) # Mediator model parameters
  sigma.resp <- par[d.resp + 1] # Error standard deviation outcome model
  sigma.expl <- par[length(par)] # Error standard deviation mediator model
  
  # Components to be used in the calculation of the gradients
  X.resp.wo2 <- X.resp[, -which(names(par)==med.name)]
  
  n <- length(outc.resp)
  
  A <- (outc.expl - X.expl%*%expl.coef)/sigma.expl
  B <- (outc.resp - X.resp%*%coef.resp)/sigma.resp
  C <- coef.resp.th2^2*sigma.expl^2 + sigma.resp^2 + 2*coef.resp.th2*Rho*sigma.expl*sigma.resp
  
  D <- (L - (X.resp.wo2%*%coef.resp.wo2 + coef.resp.th2*(X.expl%*%expl.coef)))/sqrt(C)
  
  # Calculation of the gradients
  gr.expl.coef <- colSums(as.vector(A - Rho*B)*X.expl)/((1 - Rho^2)*sigma.expl) -
    colSums(as.vector(dnorm(D)/(1 - pnorm(D)))*coef.resp.th2*X.expl/sqrt(C))
  
  
  gr.coef.resp.wo2 <- colSums(as.vector(B - Rho*A)*X.resp.wo2)/((1 - Rho^2)*sigma.resp) -
    colSums(as.vector(dnorm(D)/(1 - pnorm(D)))*X.resp.wo2)/sqrt(C)
  
  
  gr.coef.resp.th2 <- sum((B - Rho*A)*outc.expl)/((1 - Rho^2)*sigma.resp) -
       sum(dnorm(D)/(1 - pnorm(D))*(X.expl%*%expl.coef/sqrt(C) +
        (L - X.resp.wo2%*%coef.resp.wo2 - coef.resp.th2*X.expl%*%expl.coef)*(sigma.expl^2*coef.resp.th2 + 
          sigma.expl*sigma.resp*Rho)/(C^(3/2)))) 
   
  gr.sigma.b <- -n/sigma.expl + 1/(sigma.expl^2*(1 - Rho^2))*sum((outc.expl - X.expl%*%expl.coef)*(A - Rho*B)) -
    sum(dnorm(D)/(1 - pnorm(D))*((L - (X.resp.wo2%*%coef.resp.wo2 + 
       coef.resp.th2*(X.expl%*%expl.coef)))*coef.resp.th2*(coef.resp.th2*sigma.expl + sigma.resp*Rho))/(C^(3/2)))
  
  gr.sigma.t <- -n/sigma.resp + 1/(sigma.resp^2*(1 - Rho^2))*sum((outc.resp - X.resp%*%coef.resp)*(B - Rho*A)) -
    sum(dnorm(D)/(1 - pnorm(D))*((L - (X.resp.wo2%*%coef.resp.wo2 + 
     coef.resp.th2*(X.expl%*%expl.coef)))*(coef.resp.th2*Rho*sigma.expl + sigma.resp))/(C^(3/2)))
  
  gr.resp <- c(gr.coef.resp.wo2[1:2],gr.coef.resp.th2,gr.coef.resp.wo2[-c(1:2)]) 
  
  return(c(gr.resp, gr.sigma.t, gr.expl.coef, gr.sigma.b)) 
  
}

