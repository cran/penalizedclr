## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,eval = T,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("penalizedclr")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  install_github("veradjordjilovic/penalizedclr")

## ----setup--------------------------------------------------------------------
library(penalizedclr)

## ----message=FALSE, warning=FALSE---------------------------------------------
set.seed(1234)
library(tidyverse)

## -----------------------------------------------------------------------------
# two groups of predictors
p <- c(12, 40)

# percentage and number of non-null variables
p_nz <- c(0.5, 0.2)
m_nz <- round(p*p_nz, 0)

# number of different strata (case-control pairs)
K <- 100

# number of cases and controls in each stratum (not necessarily 1:1 matching,
# other designs are also allowed)
n_cases <- 1
n_ctrl <- 1


# generating covariates
X = cbind(matrix(rnorm(p[1] * K * (n_cases + n_ctrl), 0, 1), ncol = p[1]),
          matrix(rnorm(p[2] * K * (n_cases + n_ctrl), 0, 4), ncol = p[2]))

# coefficients
beta <- as.matrix(c(rnorm(m_nz[1], 4, 1),
                    rep(0, p[1] - m_nz[1]),
                    rnorm(m_nz[2], 2, 0.8),
                    rep(0, p[2] - m_nz[2])), ncol = 1)




# stratum membership
stratum <- rep(1:K, each= n_cases+n_ctrl)

# linear predictor
lin_pred <-  X %*% beta

prob_case <- exp(lin_pred) / (1 + exp(lin_pred))


# generate the response

Y <- rep(0, length(stratum))

data_sim <- as_tibble(data.frame(stratum = stratum,
                                 probability = prob_case,
                                 obs_id = 1 : length(stratum)))
data_sim_cases <- data_sim %>%
  group_by(stratum)%>%
  sample_n(n_cases, weight = probability)

Y[data_sim_cases$obs_id] <- 1

## ----results = 'hide'---------------------------------------------------------
fit1 <- penalized.clr(response = Y, penalized = X, stratum = stratum, 
                      lambda = c(1,4), p = p, standardize = TRUE)

## -----------------------------------------------------------------------------
str(fit1)
nonzero_index <- (beta != 0) * 1 #index of nonzero coefficients
table(fitted = (fit1$penalized != 0) * 1, nonzero_index)

## ----results = 'hide'---------------------------------------------------------
fit2 <- penalized.clr(response = Y, penalized = X, stratum = stratum, 
                      lambda = c(1,4),p = p, standardize = TRUE, alpha = 0.6)

## -----------------------------------------------------------------------------
stable1 <- stable.clr(response = Y, B = 50,
                      penalized = X, stratum = stratum,
                      lambda.seq = c(1, 4, 6), parallel = TRUE)

## -----------------------------------------------------------------------------
hist(stable1$Pistab, xlab = "Selection probability", ylab="Frequency")

## -----------------------------------------------------------------------------
(tab1 <- table(stable = (stable1$P>0.65) * 1, nonzero_index))
(tab1[1,1] + tab1[2,2])/(sum(tab1))


## -----------------------------------------------------------------------------
# plotting selection probabilities for true non-zero (red) and zero (blue) coefficients
s_prob_nonzero <- cut(stable1$P[nonzero_index == 1], 
                      breaks = seq(0,1, by = 0.1), ordered_result = T)
s_prob_zero <- cut(stable1$P[nonzero_index == 0], 
                   breaks = seq(0,1, by = 0.1), ordered_result = T)

barplot(table(s_prob_nonzero), col = 2)
barplot(table(s_prob_zero), add =T, col = 4 )

## ----eval= T------------------------------------------------------------------
stable2 <- penalizedclr::stable.clr.g(response = Y, p = p, 
                                      standardize = TRUE,
                                      penalized = X, stratum = stratum,
                                      lambda.list = list(c(4,1), c(1,4)))



## ----echo = FALSE, results = 'hide'-------------------------------------------
# histogram of selection probabilities
#hist(stable2$P, xlab = "Selection probability", ylab="Frequency")
# table of true status vs. selection status
(tab2 <- table(stable = (stable2$P>0.65) * 1, nonzero_index))
# classification accuracy
(tab2[1,1] + tab2[2,2])/(sum(tab2))



## ---- echo = F----------------------------------------------------------------
# plotting selection probabilities for true non-zero (red) and zero (blue) coefficients
s_prob_nonzero <- cut(stable2$P[nonzero_index == 1],
                      breaks = seq(0,1, by = 0.1), 
                      ordered_result = T)
s_prob_zero <- cut(stable2$P[nonzero_index == 0],
                   breaks = seq(0,1, by = 0.1), 
                   ordered_result = T)

barplot(table(s_prob_nonzero), col = 2)
barplot(table(s_prob_zero), add = T, col = 4 )

## ---- eval = T----------------------------------------------------------------

(pf <- default.pf(response = Y, stratum = stratum, penalized = X, 
                 nfolds = 10, alpha = 0.3,
                 standardize = TRUE, p = p, type.step1 = "comb"))




(lambda <- find.default.lambda(response = Y, penalized = X, 
                              standardize = TRUE, stratum = stratum, 
                              p = p,  pf.list  = pf))

## ----eval= T------------------------------------------------------------------
stable3 <- stable.clr.g(response = Y, p = p,  
                        standardize = TRUE,
                        penalized = X, stratum = stratum,
                        lambda.list = list(c(pf[[1]][1]*as.numeric(lambda), pf[[1]][2]*as.numeric(lambda))))

## ---- echo = F----------------------------------------------------------------
# plotting selection probabilities for true non-zero (red) and zero (blue) coefficients
s_prob_nonzero <- cut(stable3$P[nonzero_index == 1], 
                      breaks = seq(0, 1, by = 0.1), ordered_result = T)
s_prob_zero <- cut(stable3$P[nonzero_index == 0],
                   breaks = seq(0, 1, by = 0.1), ordered_result = T)

barplot(table(s_prob_nonzero), col = 2)
barplot(table(s_prob_zero), add = T, col = 4 )

## ----eval = T-----------------------------------------------------------------
alt.pf.list<- list(c(1,4), c(2,4), c(4,1), c(4,2), pf$pf)
alt.lambda <- find.default.lambda(response = Y, penalized = X, 
                              standardize = TRUE, stratum = stratum, 
                              p = p,  pf.list  = alt.pf.list)
lambda.matrix <- mapply(function (x,y) x*y, lambda, alt.pf.list)

lambda.list <- lapply(seq_len(ncol(lambda.matrix)), function(i) lambda.matrix[,i])

## -----------------------------------------------------------------------------
stable4 <- stable.clr.g(response = Y, p = p,  
                        standardize = TRUE,
                        penalized = X, stratum = stratum,
                        lambda.list = lambda.list )

## ---- echo = F----------------------------------------------------------------
# plotting selection probabilities for true non-zero (red) and zero (blue) coefficients
s_prob_nonzero <- cut(stable4$P[nonzero_index == 1], 
                      breaks = seq(0, 1, by = 0.1), ordered_result = T)
s_prob_zero <- cut(stable4$P[nonzero_index == 0],
                   breaks = seq(0, 1, by = 0.1), ordered_result = T)

barplot(table(s_prob_nonzero), col = 2)
barplot(table(s_prob_zero), add = T, col = 4 )

