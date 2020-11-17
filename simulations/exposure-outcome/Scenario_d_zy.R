#############################################
####       Simulations Scenario d        ####
####    Exposure-outcome confounding     ####
#### Binary mediator, continuous outcome ####
#############################################

# Loading required packages #

library(sensmediation)
library(mvtnorm)
library(rsimsum)

# ---------------------------- 

# Setting the values for the simulations #

n <- 500 # Sample size
S <- 2000 # Number of replicates

R <- 0.5 # True correlation between error terms of exposure and outcome model
Sigma <- cbind(c(1, R), c(R, 1)) # Covariance matrix to generate error terms for exposure and outcome model

med.coefs <- as.matrix(c(-1, 0.5, 0.5, 0.1)) # True parameters, mediator model 
out.coefs_bc <- as.matrix(c(-0.2, 0.1, 0.22, -0.15, 0.085, -0.07, -0.05, -0.04)) # True parameters, outcome model  

rho <- c(0, 0.5) # Correlation vector used for sensitivity analyses

# ---------------------------- 

# Objects to store results

nie.bc.zy.500 <- nde.bc.zy.500 <- se.nie.bc.zy.500 <- se.nde.bc.zy.500 <- matrix(nrow = S, ncol = length(rho)) 
states.bc.zy.500 <- list() # List to store the states of the random number generator, to be able to replicate results

# ----------------------------



# Simulations #

# Setting the seed:

set.seed(189690) # n = 500
# set.seed(192923) # n = 1000
# set.seed(195100) # n = 5000

tid <- proc.time()
for(i in 1:S){
  
  if(i %% 50 == 0)
    print(i)
  
  states.bc.zy.500[[i]] <- .Random.seed
  
  epsilon <- rmvnorm(n, sigma = Sigma) # Correlated error terms, exposure and outcome models
  
  X <- rnorm(n, mean = 1) # Observed confounder
  
  # Generate exposure, mediator and outcome (True models):
  Z.star <- -0.5 + 0.05*X + epsilon[, 1]
  Z <- ifelse(Z.star > 0, 1, 0)
  M.star <- med.coefs[1, ] + med.coefs[2, ]*Z + med.coefs[3, ]*X + med.coefs[4, ]*Z*X + rnorm(n)
  M <- ifelse(M.star > 0, 1, 0)
  Y <- out.coefs_bc[1, ] + out.coefs_bc[2, ]*Z + out.coefs_bc[3, ]*M + out.coefs_bc[4, ]*X + out.coefs_bc[5, ]*Z*M + out.coefs_bc[6, ]*Z*X +
    out.coefs_bc[7, ]*M*X + out.coefs_bc[8, ]*Z*M*X + epsilon[, 2]
 
  # Estimated models:
  z.model <- glm(Z ~ X, family = binomial(link = "probit"))
  m.model <- glm(M ~ Z * X, family = binomial(link = 'probit')) 
  y.model <- glm(Y ~ Z * M * X)
  
  # Estimation of effects:
  result.bc.zy.500 <- sensmediation(exp.model = z.model, med.model = m.model, out.model = y.model, Rho = rho,
                                        progress = FALSE, exp.name = "Z", med.name = "M", type = "zy")
  
  
  ###Storage of results:
  nie.bc.zy.500[i, ] <- result.bc.zy.500$NIE
  nde.bc.zy.500[i, ] <- result.bc.zy.500$NDE
  se.nie.bc.zy.500[i, ] <- result.bc.zy.500$std.errs$se.nie
  se.nde.bc.zy.500[i, ] <- result.bc.zy.500$std.errs$se.nde

  
  if(i==S)
    states.bc.zy.500[[i + 1]] <- .Random.seed
  
}
proc.time()-tid
# ---------------------------- 