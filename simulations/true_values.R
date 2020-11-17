#############################################
######## Simulations, true effects ########
#############################################

# Loading required packages #

library(mvtnorm)

# ---------------------------- 

# Setting the values for the simulations #

n <- 100000 # Sample size
S <- 10000 # Number of replicates

# True parameters
med.coefs <- as.matrix(c(-1, 0.5, 0.5, 0.1))  
out.coefs_cb <- as.matrix(c(-0.2, 0.2, 0.2, -0.15, 0.1, -0.05, -0.05, -0.05))
out.coefs_cc <- as.matrix(c(-0.2, 0.15, 0.2, -0.15, 0.1, -0.062, -0.099, -0.099))
out.coefs_bc <- as.matrix(c(-0.2, 0.1, 0.22, -0.15, 0.085, -0.07, -0.05, -0.04))
# ---------------------------- 

# Objects to store true effects

nie.True.cc <- nde.True.cc <- nie.True.bc <- nde.True.bc <- nie.True.cb <- nde.True.cb <- numeric(0)  

# ---------------------------- 


# Simulations:

# Setting the seed:

set.seed(102871)

for(i in 1:S){
  
  if(i %% 50 == 0)
    print(i)


  X <- rnorm(n, mean = 1) # Observed confounder
  
  ## Calculate effects ##
  
  # Scenario d: Binary mediator, continuous outcome
  med.ie <- pnorm(med.coefs[1, ] + med.coefs[2, ] + X*(med.coefs[3, ]  +  med.coefs[4, ]))-
    pnorm(med.coefs[1, ] + X*med.coefs[3, ])
  med.de <- pnorm(med.coefs[1, ] + X*med.coefs[3, ])
  out.ie <- out.coefs_bc[3, ]  +  out.coefs_bc[5, ]  +  X*(out.coefs_bc[7, ]  +  out.coefs_bc[8, ])
  out.de <- out.coefs_bc[5, ]  +  X*out.coefs_bc[8, ]
  nie.True.bc[i] <- mean(out.ie*med.ie)
  nde.True.bc[i] <- mean(out.coefs_bc[2, ] +  X*out.coefs_bc[6, ]  +  out.de*med.de )
  ###

  # Scenario e: Continuous mediator, binary outcome
  med.ie.tr <- med.coefs[1, ] + med.coefs[2, ] + X*(med.coefs[3, ]  +  med.coefs[4, ])
  med.ie.cont <- med.coefs[1, ] + X*med.coefs[3, ]
  med.de <- med.coefs[1, ] + X*med.coefs[3, ]

  out.ie <- out.de.tr <- out.coefs_cb[3, ]  +  out.coefs_cb[5, ] + X*(out.coefs_cb[7, ]  +  out.coefs_cb[8, ])
  out.de.cont <- out.coefs_cb[3, ]  +  X*out.coefs_cb[7, ]

  denom.ie <- sqrt(out.ie^2 + 1)

  nie <- pnorm((out.coefs_cb[1, ] + out.coefs_cb[2, ] + out.ie*med.ie.tr +  X*(out.coefs_cb[4, ] +
         out.coefs_cb[6, ]) )/denom.ie )  -
    pnorm((out.coefs_cb[1, ] + out.coefs_cb[2, ] + out.ie*med.ie.cont +  X*(out.coefs_cb[4, ] + out.coefs_cb[6, ]) )/denom.ie )

  denom.de.tr <- sqrt(out.de.tr^2 + 1)
  denom.de.cont <- sqrt(out.de.cont^2 + 1)

  nde <- pnorm((out.coefs_cb[1, ] + out.coefs_cb[2, ] + out.de.tr*med.de +  X*(out.coefs_cb[4, ] +
           out.coefs_cb[6, ]) )/denom.de.tr )  -
    pnorm((out.coefs_cb[1, ] + out.de.cont*med.de +  X*out.coefs_cb[4, ] )/denom.de.cont )

  nie.True.cb[i] <- mean(nie)
  nde.True.cb[i] <- mean(nde)
  ###

  # Scenarios a and b: Continuous mediator and outcome
  med.ie <- med.coefs[2, ] + X*med.coefs[4, ]
  med.de <- med.coefs[1, ] + X*med.coefs[3, ]

  out.ie <- out.coefs_cc[3, ] + out.coefs_cc[5, ] + X*(out.coefs_cc[7, ] + out.coefs_cc[8, ])
  out.de <- out.coefs_cc[5, ] + X*out.coefs_cc[8, ]

  nie.True.cc[i] <- mean(out.ie*med.ie)
  nde.True.cc[i] <- mean(out.coefs_cc[2, ]+ X*out.coefs_cc[6, ] + out.de*med.de )
  ###
}

# ---------------------------- 

