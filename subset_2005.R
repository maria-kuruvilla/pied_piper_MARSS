#test for hatchery and covariate effects on 2005

#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)


data <- read.csv(here("data","dungeness_subset_covariates.csv"), header = TRUE)
data$coho1_hatchery_num_interpolate <- na.approx(data$coho1_hatchery_num)

coho_2005 <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list <- list(
  U = "unequal",
  R = "diagonal and equal",
  Q = "diagonal and unequal"
)

fit <- MARSS(coho_2005, model = mod.list, method = "BFGS")

autoplot(fit)

AIC(fit)
#248

#hatchery
cov_hatchery <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  select(coho1_hatchery_num_interpolate,year,doy) %>%
  pivot_wider(names_from = year, values_from = c(coho1_hatchery_num_interpolate)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


mod.list_hatchery <- list(
  U = matrix("u"),
  R = matrix("r"),
  Q = matrix("q"),
  C = matrix("c"),
  c = cov_hatchery
)


fit_hatchery <- MARSS(coho_2005, model = mod.list_hatchery, method = "BFGS")

AIC(fit_hatchery)
#246

#temp
cov_temp <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  select(temp,year,doy) %>%
  pivot_wider(names_from = year, values_from = temp) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list_temp <- list(
  U = matrix("u"),
  R = matrix("r"),
  Q = matrix("q"),
  C = matrix("c"),
  c = cov_temp
)


fit_temp <- MARSS(coho_2005, model = mod.list_temp, method = "BFGS")

AIC(fit_temp)
#247


#flow

cov_flow <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  select(flow,year,doy) %>%
  pivot_wider(names_from = year, values_from = flow) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list_flow <- list(
  U = matrix("u"),
  R = matrix("r"),
  Q = matrix("q"),
  C = matrix("c"),
  c = cov_flow
)


fit_flow <- MARSS(coho_2005, model = mod.list_flow, method = "BFGS")

AIC(fit_flow)
#250


#photoperiod

cov_photoperiod <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  select(photoperiod,year,doy) %>%
  pivot_wider(names_from = year, values_from = photoperiod) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list_photoperiod <- list(
  U = matrix("u"),
  R = matrix("r"),
  Q = matrix("q"),
  C = matrix("c"),
  c = cov_photoperiod
)


fit_photoperiod <- MARSS(coho_2005, model = mod.list_photoperiod, method = "BFGS")

AIC(fit_photoperiod)
#247



#atu_solstice

cov_atu_solstice <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  select(atu_solstice,year,doy) %>%
  pivot_wider(names_from = year, values_from = atu_solstice) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

mod.list_atu <- list(
  U = matrix("u"),
  R = matrix("r"),
  Q = matrix("q"),
  C = matrix("c"),
  c = cov_atu_solstice
)


fit_atu <- MARSS(coho_2005, model = mod.list_atu, method = "BFGS")

AIC(fit_atu)
#247


#DLM

TT <- length(coho_2005)
m <- dim(cov_hatchery)[1] + 1

B <- diag(m)                        ## 2x2; Identity
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
Q <- matrix(list(
  "q.beta1",0,
  0, "q.beta2"), m, m)          ## 2x2; all 0 for now
# diag(Q) <- c(0, "q.beta1",  "q.beta2")   ## 2x2; diag = (q1,q2)


## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
Z[1,2,] <- cov_hatchery          ## Nx1; predictor variable
A <- matrix(0)               ## 1x1; scalar = 0
R <- matrix("r")             ## 1x1; scalar = r


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)


## fit univariate DLM
dlm_2005 <- MARSS(coho_2005, model = mod_list, inits = inits_list)

autoplot(dlm_2005)

