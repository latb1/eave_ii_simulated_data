a_k_values <- c(
  "", 
  "0", "1", 
  "00", "01", "10", "11", 
  "000", "001", "010", "011", "100", "101", "110", "111" #,
#  "0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111", "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"
)

library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
set.seed(1)
expit <- function(x) return(exp(x) / (1+exp(x)))

setwd("/Users/victorvelascopardo/eave_simulated_data/")
# setwd("/home/victor/trial_emulation_project/seaman_paper/")

demographics <- read.csv("data/multimorbidity.csv", stringsAsFactors = TRUE)
demographics$SCSIMD5 <- as.factor(demographics$SCSIMD5)
demographics$NumComorbidities <- as.factor(demographics$NumComorbidities)
# demographics$p <- demographics$n/sum(demographics$n)
demographics <- demographics[!is.na(demographics$SCSIMD5), ]

levels(demographics$AgeGroup) <- gsub("\\+", "more", gsub("-", "to", levels(demographics$AgeGroup)))

encode_binary <- function(x, name = "") {
  factors <- levels(x)
  binary_encoding <- ifelse(x == factors[2], 1, 0)
  binary_encoding <- matrix(binary_encoding, ncol = 1)
  colnames(binary_encoding) <- paste0(name, factors[2])
  return(binary_encoding)
}

encode_ordinal <- function(x, name = "") {
  factors <- levels(x)
  L <- length(factors)
  factors_matrix <- matrix(0, nrow = L, ncol = L-1)
  factors_matrix[lower.tri(factors_matrix)] <- 1
  colnames(factors_matrix) <- paste0(name, "GE", factors[2:L])
  return(factors_matrix[as.integer(x), ])
}

demographics_coded <- cbind(
  encode_binary(demographics$Sex, name = "Sex"),
  encode_ordinal(demographics$AgeGroup, name = "AgeGroup"),
  encode_ordinal(demographics$SCSIMD5, name = "SCSIMD5"),
  encode_ordinal(demographics$NumComorbidities, name = "NumComorbidities"),
  demographics$n
)
colnames(demographics_coded)[ncol(demographics_coded)] <- "n"

# browser()
# Set up N (sample size), K (number of visits)
m <- 1e5 # Sample size for estimating the CDF
K <- 2    # Number of visits (minus one)

# Parameters for causal quantity of interest
betas <- list()
betas[[1]] <- c(-2.5, -0.5)
betas[[2]] <- c(-2.5, -0.5, -0.5)
betas[[3]] <- c(-2.5, 0, -0.5, -0.5)

# Parameters for allocation of treatment to individuals
gammas <- list()
gammas[[1]] <- c(
  -1.5,
  0.05,
  rep(0.1, nlevels(demographics$AgeGroup)-1),
  rep(0.05, nlevels(demographics$SCSIMD5)-1),
  rep(0.1, nlevels(demographics$NumComorbidities)-1)
)
gammas[[2]] <- c(
  -1.5,
  0.05,
  rep(0.1, nlevels(demographics$AgeGroup)-1),
  rep(0.05, nlevels(demographics$SCSIMD5)-1),
  rep(0.1, nlevels(demographics$NumComorbidities)-1),
  -0.5
)
gammas[[3]] <- c(
  -1.5,
  0.05,
  rep(0.1, nlevels(demographics$AgeGroup)-1),
  rep(0.05, nlevels(demographics$SCSIMD5)-1),
  rep(0.1, nlevels(demographics$NumComorbidities)-1),
  0, -0.5
)
# gammas[[4]] <- c(gammas[[1]], c(0, -0.5, -0.5))

# Parameters for confounding mechanism
thetas <- c(
  0.05,
  rep(0.1, nlevels(demographics$AgeGroup)-1),
  rep(0.05, nlevels(demographics$SCSIMD5)-1),
  rep(0.1, nlevels(demographics$NumComorbidities)-1)
)
rho <- rep(-0.1, K+1)

fhat <- list()

for (tmt_allocation in a_k_values) {
  # Step 1. For j=1,...,m, simulate B_j ~ p(B_j). Set k = 0.
  b <- demographics[sample(1:nrow(demographics), size = m, replace = TRUE, prob = demographics$n), -ncol(demographics)]
  rownames(b) <- NULL
  b_coded <- cbind(
    encode_binary(b$Sex, name = "Sex"),
    encode_ordinal(b$AgeGroup, name = "AgeGroup"),
    encode_ordinal(b$SCSIMD5, name = "SCSIMD5"),
    encode_ordinal(b$NumComorbidities, name = "NumComorbidities")
  )
  for (k in 0:nchar(tmt_allocation)) {
    cat("Tmt allocation = ", tmt_allocation, " k = ", k, "\n")
    # (Skip) Step 2. Sample time-varying covariates
    # Step 3. Compute H_k, \hat{F}_H.
    h <- as.vector(b_coded %*% thetas)
    if (tmt_allocation == "") {
      fhat[["BASE"]] <- ecdf(h)
    } else if (k == nchar(tmt_allocation)) {
      fhat[[tmt_allocation]] <- ecdf(h)
    }
    
    # Step 4. If k = K: Stop
    if (k == nchar(tmt_allocation)) break
    
    # Step 5. Compute R_k, ranking of H_k
    r <- rank(h)
    
    # Step 6. Computer U_H
    u_h <- (r - runif(m))/m
    
    # Step 7. Compute Z_H
    z_h <- qnorm(u_h)
    
    # Step 8. Sample Z_Y, calculate U_Y
    z_y <- rnorm(m, rho[k+1]*z_h, sqrt(1-rho[k+1]^2))
    u_y <- pnorm(z_y)
    
    # Step 9. Compute Y_{k+1}
    failed <- u_y < expit(sum(as.vector(c(1, as.numeric(strsplit(substr(tmt_allocation, 1, k+1), "")[[1]]))) * betas[[k+1]]))
    
    # Step 10. Replace failed individuals with individuals that have survived
    random_replacements <- sample(which(!failed), size = sum(failed), replace = TRUE)
    b[which(failed), ] <- b[random_replacements, ]
    b_coded[which(failed), ] <- b_coded[random_replacements, ]
    
    # Step 11. Set k = k + 1 and return to step 2
    
  }
}

