#############################################
######     Simulations Scenario a      ######
######   Exposure-outcome confounding  ######
###### Continuous mediator and outcome ######
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
Sigma <- cbind(c(1, R), c(R, 1)) # Covariance matrix to generate error terms for exposure and outcome models

med.coefs <- as.matrix(c(-1, 0.5, 0.5, 0.1)) # True parameters, mediator model  
out.coefs_cc <- as.matrix(c(-0.2, 0.15, 0.2, -0.15, 0.1, -0.062, -0.099, -0.099)) # True parameters, outcome model   

rho <- c(0, 0.5) # Correlation vector used for sensitivity analyses

# ---------------------------- 

# Objects to store results

nie.cc.zy.500 <- nde.cc.zy.500 <- se.nie.cc.zy.500 <- se.nde.cc.zy.500 <- matrix(nrow = S, ncol = length(rho)) 
states.cc.zy.500 <- list() # List to store the states of the random number generator, to be able to replicate results

# ----------------------------



# Simulations #

# Setting the seed:

set.seed(180667) # n = 500
# set.seed(183622) # n = 1000
# set.seed(184765) # n = 5000

tid <- proc.time()
for(i in 1:S){
  
  if(i %% 50 == 0)
    print(i)
  
  states.cc.zy.500[[i]] <- .Random.seed
  
  epsilon <- rmvnorm(n, sigma = Sigma) # Correlated error terms, exposure and outcome models
  
  X <- rnorm(n, mean = 1) # Observed confounder
  
  # Generate exposure, mediator and outcome (True models):
  Z.star <- -0.5 + 0.05*X + epsilon[, 1] 
  Z <- ifelse(Z.star > 0, 1, 0)
  M <- med.coefs[1, ] + med.coefs[2, ]*Z + med.coefs[3, ]*X + med.coefs[4, ]*Z*X + rnorm(n)
  Y <- out.coefs_cc[1, ] + out.coefs_cc[2, ]*Z + out.coefs_cc[3, ]*M + out.coefs_cc[4, ]*X + out.coefs_cc[5, ]*Z*M + out.coefs_cc[6, ]*Z*X +
    out.coefs_cc[7, ]*M*X + out.coefs_cc[8, ]*Z*M*X + epsilon[, 2] 
  
  # Estimated models:
  z.model <- glm(Z ~ X, family = binomial(link = "probit"))
  m.model <- glm(M ~ Z * X) 
  y.model <- glm(Y ~ Z * M * X)
  
  # Estimation of effects:
  result.cc.zy.500 <- sensmediation(exp.model = z.model, med.model = m.model, out.model = y.model, Rho = rho,
                                        progress = FALSE, exp.name = "Z", med.name = "M", type = "zy")
  
  
  ###Storage of results:
  nie.cc.zy.500[i, ] <- result.cc.zy.500$NIE
  nde.cc.zy.500[i, ] <- result.cc.zy.500$NDE
  se.nie.cc.zy.500[i, ] <- result.cc.zy.500$std.errs$se.nie
  se.nde.cc.zy.500[i, ] <- result.cc.zy.500$std.errs$se.nde
  
  
  if(i==S)
    states.cc.zy.500[[i + 1]] <- .Random.seed
  
}
proc.time()-tid
# ---------------------------- 
