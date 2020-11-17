#############################################
######     Simulations Scenario b      ######
######   Mediator-outcome confounding  ######
###### Continuous mediator and outcome ######
######       Censored outcome          ######
#############################################

# Loading required packages #

library(sensmediation)
library(mvtnorm)
library(rsimsum)
library(AER)

# ---------------------------- 

# Setting the values for the simulations #

n <- 500 # Sample size
S <- 2000 # Number of replicates

R <- 0.5 # True correlation between error terms of mediator and outcome model
Sigma <- cbind(c(1, R), c(R, 1)) # Covariance matrix to generate error terms for mediator and outcome models

med.coefs <- as.matrix(c(-1, 0.5, 0.5, 0.1)) # True parameters, mediator model  
out.coefs_cc <- as.matrix(c(-0.2, 0.15, 0.2, -0.15, 0.1, -0.062, -0.099, -0.099)) # True parameters, outcome model   

rho <- c(0, 0.5) # Correlation vector used for sensitivity analyses

# ---------------------------- 

# Objects to store results

nie.cens.my.500 <- nde.cens.my.500 <- se.nie.cens.my.500 <- se.nde.cens.my.500 <- matrix(nrow = S, ncol = length(rho)) 
states.cens.my.500 <- list() # List to store the states of the random number generator, to be able to replicate results

# ----------------------------



# Simulations #

# Setting the seed:

set.seed(64553) # n = 500
# set.seed(66569) # n = 1000
# set.seed(68581) # n = 5000

for(i in 1:S){
  
  if(i %% 50 == 0)
    print(i)
  
  states.cens.my.500[[i]] <- .Random.seed
  
  epsilon <- rmvnorm(n, sigma = Sigma) # Correlated error terms, mediator and outcome models
  
  X <- rnorm(n, mean = 1) # Observed confounder
  
  # Generate exposure, mediator and outcome (True models):
  Z.star <- -0.5 + 0.05*X + rnorm(n)
  Z <- ifelse(Z.star > 0, 1, 0)
  M <- med.coefs[1, ] + med.coefs[2, ]*Z + med.coefs[3, ]*X + med.coefs[4, ]*Z*X + epsilon[, 1]
  Y <- out.coefs_cc[1, ] + out.coefs_cc[2, ]*Z + out.coefs_cc[3, ]*M + out.coefs_cc[4, ]*X + out.coefs_cc[5, ]*Z*M + out.coefs_cc[6, ]*Z*X +
    out.coefs_cc[7, ]*M*X + out.coefs_cc[8, ]*Z*M*X + epsilon[, 2] 
  Y <- ifelse(Y <= -1.3, -1.3, Y) # Censoring Y (around 20% censoring) 


  # Estimated models:
  m.model <- glm(M ~ Z * X)
  y.model <- AER::tobit(Y ~ Z * M * X, left = -1.3, right = Inf, dist = "gaussian") 

  # Estimation of effects:
  MLcoefs <- ML.cens(m.model, y.model, Rho = rho, L = -1.3, progress = FALSE) 
  result.cens.my.500 <- calc.effects(MLcoefs, exp.name = "Z", med.name = "M")

  ###Storage of results:
  nie.cens.my.500[i, ] <- result.cens.my.500$effects$NIE
  nde.cens.my.500[i, ] <- result.cens.my.500$effects$NDE
  se.nie.cens.my.500[i, ] <- result.cens.my.500$std.errs$se.nie
  se.nde.cens.my.500[i, ] <- result.cens.my.500$std.errs$se.nde


  if(i==S)
    states.cens.my.500[[i + 1]] <- .Random.seed
  
}

# ---------------------------- 