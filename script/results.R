library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

setwd("/Users/victorvelascopardo/eave_ii_simulated_data/")
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
  
  b <- read_csv(paste0("simulations/", l, "/confounders.csv")) |>
    mutate_all(factor) |>
    as.data.frame()
  b_coded <- cbind(
    encode_binary(b$Sex, name = "Sex"),
    encode_ordinal(b$AgeGroup, name = "AgeGroup"),
    encode_ordinal(b$SCSIMD5, name = "SCSIMD5"),
    encode_ordinal(b$NumComorbidities, name = "NumComorbidities")
  ) 
  a <- read_csv(paste0("simulations/", l, "/treatment.csv")) |> as.data.frame() |> as.matrix()
  y <- read_csv(paste0("simulations/", l, "/outcome.csv")) |> as.data.frame() |> as.matrix()
  
  failure.times <- rowSums(y)
  for (k in 1:(K+1)) {
    model.outcome.unweighted[[k]] <- glm(
      as.factor(y[,k+1]==0) ~ a[,1:k], 
      family = binomial, 
      subset = which(failure.times>=k)
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
      fit.numer <- glm(a[, k] ~ a[,(1:(k-1))], family = binomial, subset = which(failure.times>=k))
      fit.denom <- glm(a[, k] ~ a[,(1:(k-1))] + b_coded, family = binomial, subset = which(failure.times>=k))
    }
    
    phat.denom <- a[failure.times>=k, k] * fitted(fit.denom) + (1-a[failure.times>=k, k]) * (1- fitted(fit.denom))
    phat.numer <- a[failure.times>=k, k] * fitted(fit.numer) + (1-a[failure.times>=k, k]) * (1- fitted(fit.numer))
    if (k == 1) {
      w[, k] <- phat.numer / phat.denom
    } else {
      w[failure.times>=k, k] <- w[failure.times>=k,k-1] * phat.numer / phat.denom
    }
    model.outcome.weighted[[k]] <- glm(
      as.factor(y[,k+1]==0) ~ a[,1:k], 
      family = binomial, 
      weights = w[, k], 
      subset = which(failure.times>=k)
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
  write_csv(as.data.frame(param.estimates), paste0("output/", l, "/results.csv"))
}
