#############################################
####       Simulations Scenario e        ####
####    Exposure-outcome confounding     ####
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

R <- 0.5 # True correlation between error terms of exposure and outcome models
Sigma <- cbind(c(1, R), c(R, 1)) # Covariance matrix to generate error terms for exposure and outcome models

med.coefs <- as.matrix(c(-1, 0.5, 0.5, 0.1)) # True parameters, mediator model  
out.coefs_cb <- as.matrix(c(-0.2, 0.2, 0.2, -0.15, 0.1, -0.05, -0.05, -0.05)) # True parameters, outcome model   

rho <- c(0, 0.5) # Correlation vector used for sensitivity analyses

# ---------------------------- 

# Objects to store results

nie.cb.zy.500 <- nde.cb.zy.500 <- se.nie.cb.zy.500 <- se.nde.cb.zy.500 <- matrix(nrow = S, ncol = length(rho)) 
states.cb.zy.500 <- list() # List to store the states of the random number generator, to be able to replicate results

# ---------------------------- 


# Simulations:

# Setting the seed:

set.seed(200949) # n = 500
# set.seed(203351) # n = 1000
# set.seed(205550) # n = 5000

tid <- proc.time()
for(i in 1:S){
  
  if(i %% 50 == 0)
    print(i)
  
  states.cb.zy.500[[i]] <- .Random.seed
  
  epsilon <- rmvnorm(n, sigma = Sigma) # Correlated error terms, exposure and outcome models
  
  X <- rnorm(n, mean = 1) # Observed confounder
  
  # Generate exposure, mediator and outcome (True models):
  Z.star <- -0.5 + 0.05*X + epsilon[, 1]
  Z <- ifelse(Z.star > 0, 1, 0)
  M <- med.coefs[1, ] + med.coefs[2, ]*Z + med.coefs[3, ]*X + med.coefs[4, ]*Z*X + rnorm(n)
  Y.star <- out.coefs_cb[1, ] + out.coefs_cb[2, ]*Z + out.coefs_cb[3, ]*M + out.coefs_cb[4, ]*X + out.coefs_cb[5, ]*Z*M + out.coefs_cb[6, ]*Z*X +
    out.coefs_cb[7, ]*M*X + out.coefs_cb[8, ]*Z*M*X + epsilon[, 2]
  Y <- ifelse(Y.star > 0, 1, 0)
  
  # Estimated models:
  z.model <- glm(Z ~ X, family = binomial(link = "probit"))
  m.model <- glm(M ~ Z * X) 
  y.model <- glm(Y ~ Z * M * X, family = binomial(link = 'probit'))
  
  # Estimation of effects:
  result.cb.zy.500 <- sensmediation(exp.model = z.model, med.model = m.model, out.model = y.model, Rho = rho,
                                    progress = FALSE, exp.name = "Z", med.name = "M", type = "zy")
  
 
  ###Storage of results:
  nie.cb.zy.500[i, ] <- result.cb.zy.500$NIE
  nde.cb.zy.500[i, ] <- result.cb.zy.500$NDE
  se.nie.cb.zy.500[i, ] <- result.cb.zy.500$std.errs$se.nie
  se.nde.cb.zy.500[i, ] <- result.cb.zy.500$std.errs$se.nde
  
  
  if(i==S)
    states.cb.zy.500[[i + 1]] <- .Random.seed
  
}
proc.time()-tid
# ---------------------------- 