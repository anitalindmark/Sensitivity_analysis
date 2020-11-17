#############################################
####       Simulations Scenario e        ####
####     Mediator-outcome confounding    ####
#### Continuous mediator, binary outcome ####
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
out.coefs_cb <- as.matrix(c(-0.2, 0.2, 0.2, -0.15, 0.1, -0.05, -0.05, -0.05)) # True parameters, outcome model   

rho <- c(0, 0.5) # Correlation vector used for sensitivity analyses

# ---------------------------- 

# Objects to store results

nie.cb.my.500 <- nde.cb.my.500 <- se.nie.cb.my.500 <- se.nde.cb.my.500 <- matrix(nrow = S, ncol = length(rho)) 
states.cb.my.500 <- list() # List to store the states of the random number generator, to be able to replicate results

# ---------------------------- 


# Simulations:

# Setting the seed:

set.seed(48491) # n = 500
# set.seed(50497) # n = 1000
# set.seed(52501) # n = 5000

for(i in 1:S){
  
  if(i %% 50 == 0)
    print(i)
  
  states.cb.my.500[[i]] <- .Random.seed
  
  epsilon <- rmvnorm(n, sigma = Sigma) # Correlated error terms, mediator and outcome models
  
  X <- rnorm(n, mean = 1) # Observed confounder
  
  # Generate exposure, mediator and outcome (True models):
  Z.star <- -0.5 + 0.05*X + rnorm(n)
  Z <- ifelse(Z.star > 0, 1, 0)
  M <- med.coefs[1, ] + med.coefs[2, ]*Z + med.coefs[3, ]*X + med.coefs[4, ]*Z*X + epsilon[, 1]
  Y.star <- out.coefs_cb[1, ] + out.coefs_cb[2, ]*Z + out.coefs_cb[3, ]*M + out.coefs_cb[4, ]*X + out.coefs_cb[5, ]*Z*M + out.coefs_cb[6, ]*Z*X +
    out.coefs_cb[7, ]*M*X + out.coefs_cb[8, ]*Z*M*X + epsilon[, 2]
  Y <- ifelse(Y.star > 0, 1, 0)
  
  # Estimated models:
  m.model <- glm(M ~ Z * X) 
  y.model <- glm(Y ~ Z * M * X, family = binomial(link = 'probit'))
  
  # Estimation of effects:
  result.cb.my.500 <- sensmediation(med.model = m.model, out.model = y.model, Rho = rho,
                                      progress = FALSE, exp.name = "Z", med.name = "M")
  
 
  ###Storage of results:
  nie.cb.my.500[i, ] <- result.cb.my.500$NIE
  nde.cb.my.500[i, ] <- result.cb.my.500$NDE
  se.nie.cb.my.500[i, ] <- result.cb.my.500$std.errs$se.nie
  se.nde.cb.my.500[i, ] <- result.cb.my.500$std.errs$se.nde
  
  
  if(i==S)
    states.cb.my.500[[i + 1]] <- .Random.seed
  
}

# ---------------------------- 
