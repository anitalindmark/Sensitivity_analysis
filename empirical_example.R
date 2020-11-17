########################################
######     Empirical example      ######
######   using the UPB-data from  ######
######    the medflex package     ######
########################################

# Loading the data #
library(medflex)
data("UPBdata")



# Specification of mediator, outcome and exposure models #
medmod <- glm(negaff ~ attbin * gender + educ + age, data = UPBdata, family = gaussian)

outmod <- glm(UPB ~ (attbin * negaff)* gender + educ + age, data = UPBdata, family = binomial(link = "probit"))

expmod <- glm(attbin ~ gender + educ + age, data = UPBdata, family = binomial(link = "probit"))


## Mediation and sensitivity analyses ##
library(sensmediation)

rhos <- seq(-0.9,0.9,0.1) # Sensitivity parameter vector

# Marginal effects objects #
effects_MY <- sensmediation(med.model = medmod, out.model = outmod,
                               exp.name = "attbin", med.name = "negaff", type = "my",
                               Rho = rhos)

effects_ZM <- sensmediation(med.model = medmod, out.model = outmod,
                            exp.model = expmod, exp.name = "attbin", med.name = "negaff",
                            type = "zm", Rho = rhos)


effects_ZY <- sensmediation(med.model = medmod, out.model = outmod,
                            exp.model = expmod, exp.name = "attbin", med.name = "negaff",
                            type = "zy", Rho = rhos)

# Summary of results and plots are obtained by calling summary() and plot() on an object

# # Conditional effects objects #
effects_MY_Male <- more.effects(effects.MY, covariates = list(gender = "M"))
effects_ZM_Male <- more.effects(effects.ZM, covariates = list(gender = "M"))
effects_ZY_Male <- more.effects(effects.ZY, covariates = list(gender = "M"))

effects_MY_Female <- more.effects(effects.MY, covariates = list(gender = "F"))
effects_ZM_Female <- more.effects(effects.ZM, covariates = list(gender = "F"))
effects_ZY_Female <- more.effects(effects.ZY, covariates = list(gender = "F"))

