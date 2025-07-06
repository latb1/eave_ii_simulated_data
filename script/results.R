library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

setwd("C:/Users/the_o/Desktop/Dissertation/SynthesisationRepo/eave_ii_simulated_data")
K <- 2 # Number of visits minus one

# Parameter determining the strength of confounding (five scenarios)
rhos <- c(-0.1, -0.3, -0.5, -0.7, -0.9)

source("script/aux_functions.R")
# Parameters for causal quantity of interest
betas <- list()
betas[[1]] <- c(-3.5, -0.5)
betas[[2]] <- c(-3.5, -0.5, -0.5)
betas[[3]] <- c(-3.5, 0, -0.5, -0.5)

for (l in 1:length(rhos)) {
  model.outcome.unweighted <- list()
  rho <- rhos[l]
  
  # Extract data from the df passed from main.R
  n_confounders <- 4  # Sex, AgeGroup, SCSIMD5, NumComorbidities
  n_treatments <- 3   # A0, A1, A2

  b <- df[, 1:n_confounders] |>
    as.data.frame()
  
  # Assign proper column names since they were lost in Python conversion
  colnames(b) <- c("Sex", "AgeGroup", "SCSIMD5", "NumComorbidities")
  
  # Convert to factors
  b <- b |> mutate_all(factor)
  
  a <- df[, (n_confounders + 1):(n_confounders + n_treatments)] |> 
    as.data.frame() |> 
    as.matrix()
  y <- df[, (n_confounders + n_treatments + 1):ncol(df)] |> 
    as.data.frame() |> 
    as.matrix()
    
  # Now b_coded will work because column names exist
  b_coded <- cbind(
    encode_binary(b$Sex, name = "Sex"),
    encode_ordinal(b$AgeGroup, name = "AgeGroup"),
    encode_ordinal(b$SCSIMD5, name = "SCSIMD5"),
    encode_ordinal(b$NumComorbidities, name = "NumComorbidities")
  )
  
  failure_times <- rowSums(y)
  for (k in 1:(K+1)) {
    model.outcome.unweighted[[k]] <- glm(
      as.factor(y[,k+1]==0) ~ a[,1:k], 
      family = binomial, 
      subset = which(failure_times>=k)
    ) 
  }
  
  N <- nrow(b)
  w <- matrix(1, N, K+1)
  model.outcome.weighted <- list()
  for (k in 1:(K+1)) {
    if (k == 1) {
      fit.numer <- glm(a[, 1] ~ 1, family = binomial)
      fit.denom <- glm(a[, 1] ~ b_coded, family = binomial)
    } else {
      fit.numer <- glm(a[, k] ~ a[,(1:(k-1))], family = binomial, subset = which(failure_times>=k))
      fit.denom <- glm(a[, k] ~ a[,(1:(k-1))] + b_coded, family = binomial, subset = which(failure_times>=k))
    }
    
    phat.denom <- a[failure_times>=k, k] * fitted(fit.denom) + (1-a[failure_times>=k, k]) * (1- fitted(fit.denom))
    phat.numer <- a[failure_times>=k, k] * fitted(fit.numer) + (1-a[failure_times>=k, k]) * (1- fitted(fit.numer))
    if (k == 1) {
      w[, k] <- phat.numer / phat.denom
    } else {
      w[failure_times>=k, k] <- w[failure_times>=k,k-1] * phat.numer / phat.denom
    }
    model.outcome.weighted[[k]] <- glm(
      as.factor(y[,k+1]==0) ~ a[,1:k], 
      family = binomial, 
      weights = w[, k], 
      subset = which(failure_times>=k)
    )
  }
  
  param.estimates <- do.call(rbind, lapply(1:3, 
    function(i) cbind(
      betas[[i]],
      coefficients(summary(model.outcome.unweighted[[i]]))[,1:2],
      coefficients(summary(model.outcome.weighted[[i]]))[,1:2]
    ) |> round(4)
  ))
  
  colnames(param.estimates)[1] <- "True value"
  # write_csv(as.data.frame(param.estimates), paste0("output/", l, "/results.csv"))

  cat("\n=== RESULTS FOR RHO =", rhos[l], "===\n")
  colnames(param.estimates) <- c("True value", "Unweighted Estimate", "SD", "Weighted Estimate", "SD")
  rownames(param.estimates) <- paste0("beta_", c("00", "01", "10", "11", "12", "20", "21", "22", "23"))
  print(param.estimates)
  cat("\n")
}
