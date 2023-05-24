### Goal - try DLM for salmon smolt migration

#packages
library(here)
library(MARSS)
library(tidyverse)
library(ggplot2)
library(zoo)


#####
# survival

## load the data
data(SalmonSurvCUI, package = "MARSS")
## get time indices
years <- SalmonSurvCUI[, 1]
## number of years of data
TT <- length(years)
## get response variable: logit(survival)
dat <- matrix(SalmonSurvCUI[,2], nrow = 1)
#####

####### migration

data <- read.csv(here("data","dungeness_all_days.csv"), header = TRUE)

days = data$doy[data$year == 2005]

#filtering only 2005-2009 data for coho wild

coho <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

TT <- length(coho)


#####

# survivial

## get predictor variable
CUI <- SalmonSurvCUI[,3]
## z-score the CUI
CUI_z <- matrix((CUI - mean(CUI)) / sqrt(var(CUI)), nrow = 1)
## number of regr params (slope + intercept)
m <- dim(CUI_z)[1] + 1



######


### migration
data$coho1_hatchery_num_interpolate <- na.approx(data$coho1_hatchery_num)

cov <- arrange(data,doy) %>%
  filter(year == 2005 & doy > 100 & doy <=200) %>%
  select(coho1_hatchery_num_interpolate,year,doy, temp) %>%
  pivot_wider(names_from = year, values_from = c(coho1_hatchery_num_interpolate,temp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


m_cov <- dim(cov)[1] + 1


## for process eqn
B <- diag(m_cov)                        ## 2x2; Identity
U <- matrix(0, nrow = m_cov, ncol = 1)  ## 2x1; both elements = 0
Q <- matrix(list(
  0,0,0,
  0, "q.beta1",0,
  0,0,"q.beta2"), m_cov, m_cov)          ## 2x2; all 0 for now
# diag(Q) <- c(0, "q.beta1",  "q.beta2")   ## 2x2; diag = (q1,q2)



## for observation eqn
Z <- array(NA, c(1, m_cov, TT))  ## NxMxT; empty for now
Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
Z[1,2:3,] <- cov            ## Nx1; predictor variable
A <- matrix(0)               ## 1x1; scalar = 0
R <- matrix("r")             ## 1x1; scalar = r



## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m_cov))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)


## fit univariate DLM
dlm_new <- MARSS(coho, model = mod_list, inits = inits_list)

plot(dlm_new$states[1,], type = "l", ylim = c(0,5))
lines(dlm_new$states[1,] + dlm_new$states.se[1,], lty = 2)
lines(dlm_new$states[1,] - dlm_new$states.se[1,], lty = 2)
abline(h=0)

plot(dlm_new$states[2,], type = "l", ylim = c(-1,1))
lines(dlm_new$states[2,] + dlm_new$states.se[2,], lty = 2)
lines(dlm_new$states[2,] - dlm_new$states.se[2,], lty = 2)
abline(h=0)

plot(dlm_new$states[3,], type = "l", ylim = c(-1,1))

lines(dlm_new$states[3,] + dlm_new$states.se[3,], lty = 2)
lines(dlm_new$states[3,] - dlm_new$states.se[3,], lty = 2)

abline(h=0)


######
## trying with multiple years

coho <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2008 & doy > 100 & doy <=200) %>%
  mutate(log.value = ifelse(coho1_wild_num == 0, NA, log(coho1_wild_num))) %>%
  select(log.value,year,doy) %>%
  pivot_wider(names_from = year, values_from = log.value) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()

TT <- length(coho[1,])

data$coho1_hatchery_num_interpolate <- na.approx(data$coho1_hatchery_num)

cov <- arrange(data,doy) %>%
  filter(year >= 2005 & year < 2008 & doy > 100 & doy <=200) %>%
  select(coho1_hatchery_num_interpolate,year,doy, temp) %>%
  pivot_wider(names_from = year, values_from = c(coho1_hatchery_num_interpolate,temp)) %>%
  column_to_rownames(var = "doy") %>%
  as.matrix() %>%
  t()


m_cov <- dim(cov)[1]/dim(coho)[1] + 1


## for process eqn
B <- diag(m_cov)                        ## 2x2; Identity
U <- matrix(0, nrow = m_cov, ncol = 1)  ## 2x1; both elements = 0
# Q <- array(NA, c(dim(coho)[1], m_cov, m_cov))
# Q[1,1,] <- rep(0,m_cov)
# Q[2,1,] <- rep(0,m_cov)
# Q[3,1,] <- rep(0,m_cov)
# 
# Q[1,2,] <- matrix(list(0,"q1",0, 0, 0, "q2"))
# Q[2,2,] <- rep(0,m_cov)
# Q[3,2,] <- rep(0,m_cov)


# Q <- matrix(list(
#   0,0,0,0,0,0,0,
#   0,"q.h2005",0,0,0,0,0,
#   0,0,"q.h2006",0,0,0,0,
#   0,0,0,"q.h2007",0,0,0,
#   0,0,0,0,"q.t2005",0,0,
#   0,0,0,0,0,"q.t2006",0,
#   0,0,0,0,0,0,"q.t2007"), m_cov, m_cov)          ## 2x2; all 0 for now
# # diag(Q) <- c(0, "q.beta1",  "q.beta2")   ## 2x2; diag = (q1,q2)



## for observation eqn
Z <- array(NA, c(3, 3, TT))  ## NxMxT; empty for now
Z[1,1,] <- rep(1, TT)        ## Nx1; 1's for intercept
Z[1,2,] <- cov[1,] ## Nx1; predictor variable
Z[1,3,] <- cov[4,] ## Nx1; predictor variable
Z[2,1,] <-  rep(1, TT) 
Z[2,2,] <- cov[2,] 
Z[2,3,] <- cov[5,] 
Z[3,1,] <- rep(1, TT) 
Z[3,2,] <- cov[3,] 
Z[3,3,] <- cov[6,] 


A <- matrix(list(0,0,0))  
R <- matrix(0, nrow = m_cov, ncol = m_cov) ## 1x1; scalar = 0
diag(R) <- rep("r", 3)             ## 1x1; scalar = r



## only need starting values for regr parameters
inits_list <-  list(x0 = matrix(c(0, 0, 0), nrow = m_cov))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = "diagonal and unequal", Z = Z, A = A, R = R)


## fit multiivariate DLM
dlm_multi <- MARSS(coho, model = mod_list, inits = inits_list)
