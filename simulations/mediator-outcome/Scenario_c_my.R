#############################################
######     Simulations Scenario c      ######
######   Mediator-outcome confounding  ######
###### Continuous mediator and outcome ######
######        Truncated outcome        ######
#############################################

# Loading required packages #

library(sensmediation)
library(mvtnorm)
library(rsimsum)

# ---------------------------- 

# Setting the values for the simulations #

n <- 500 # Sample size
S <- 2000 # Number of replicates

R <- 0.5 # True correlation between error terms of mediator and outcome model
Sigma <- cbind(c(1, R), c(R, 1)) # Covariance matrix to generate error terms for mediator and outcome models

med.coefs <- as.matrix(c(-1, 0.5, 0.5, 0.1)) # True parameters, mediator model
out.coefs_cc_main <- as.matrix(c(-1, 0.038, 0.082, -0.3)) # True parameters, outcome model

rho <- c(0, 0.5) # Correlation vector used for sensitivity analyses

# ---------------------------- 

# Objects to store results

nie.trunc.my.500 <- nde.trunc.my.500 <- se.nie.trunc.my.500 <- se.nde.trunc.my.500 <- matrix(nrow = S, ncol = length(rho)) 
states.trunc.my.500 <- list() # List to store the states of the random number generator, to be able to replicate results

# ----------------------------



# Simulations #

# Setting the seed:

set.seed(32401) # n = 500
# set.seed(34403) # n = 1000
# set.seed(36433) # n = 5000

tid <- proc.time()
for(i in 1:S){
  
  if(i %% 50 == 0)
    print(i)

  states.trunc.my.500[[i]] <- .Random.seed
  
  epsilon <- rmvnorm(n, sigma = Sigma) # Correlated error terms, mediator and outcome models
  
  X <- rnorm(n, mean = 1) # Observed confounder
  
  # Generate exposure, mediator and outcome (True models):
  Z.star <- -0.5 + 0.05*X + rnorm(n)
  Z <- ifelse(Z.star > 0, 1, 0)
  M <- med.coefs[1, ] + med.coefs[2, ]*Z + med.coefs[3, ]*X + epsilon[, 1]
  Y <- out.coefs_cc_main[1, ] + out.coefs_cc_main[2, ]*Z + out.coefs_cc_main[3, ]*M + out.coefs_cc_main[4, ]*X + epsilon[, 2] 
  datat <- data.frame(Y, M, Z, X)
  datat.trunc <- subset(datat, Y > -2.2) # truncating the data (around 20% truncation)

  # # Estimated models:
  m.model <- glm(M ~ Z + X, data = datat.trunc)
  y.model <- glm(Y ~ Z + M + X, data = datat.trunc)


  # Estimation of effects:
  MLcoefs <- ML.trunc(m.model, y.model, Rho = rho, L = -2.2, progress = FALSE, med.name = "M")
  result.trunc.my.500 <- calc.effects(MLcoefs, exp.name = "Z", med.name = "M")

  # Storage of results
  nie.trunc.my.500[i, ] <- result.trunc.my.500$effects$NIE
  nde.trunc.my.500[i, ] <- result.trunc.my.500$effects$NDE
  se.nie.trunc.my.500[i, ] <- result.trunc.my.500$std.errs$se.nie
  se.nde.trunc.my.500[i, ] <- result.trunc.my.500$std.errs$se.nde

  if(i==S)
    states.trunc.my.500[[i+1]] <- .Random.seed
}
proc.time()-tid
# ---------------------------- 

