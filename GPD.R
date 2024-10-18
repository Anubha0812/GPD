# ---------------------------------------------------------------------------------------------------
# GPD Parameter Estimation and Hypothesis Testing between Pre and Post Announcement Periods
# ---------------------------------------------------------------------------------------------------
# This script compares the Gibbs Distribution parameter estimates of the persistence diagrams for transaction data
# in the pre-announcement and post-announcement periods for multiple agents. The output includes:
# 1. Parameter estimates for the GPD model in both periods.
# 2. P-values from hypothesis testing between parameters from the two periods.
# 3. A combined matrix of p-values, parameter estimates, and variances.
#
# INPUT:
# - input_pre: list of transaction data for agents in the pre-announcement period.
# - input_post: list of transaction data for agents in the post-announcement period.
#
# OUTPUT:
# - Parameters for GPD models in both periods.
# - P-values for hypothesis testing.
# - Combined matrix of p-values, parameter estimates, and variances.
# ---------------------------------------------------------------------------------------------------

# Libraries
library(foreach)
library(doParallel)
library(parallel)
library(TDA)
library(pracma)
library(dbscan)
library(ucminf)
library(mvtnorm)
library(ks)

# Setup
K <- 3  # Number of GPD parameters excluding alpha
Z <- 40  # Minimum number of transactions required to estimate the GPD model
eta <- 0.1  # Bandwidth parameter for KDE

input_pre <- list()  # Pre-announcement transaction data for agents
input_post <- list()  # Post-announcement transaction data for agents

# Storage for results
DiagGrid1 <- list()  # Persistence diagrams for pre-announcement period
DiagGrid2 <- list()  # Persistence diagrams for post-announcement period
pval <- list()  # P-values for hypothesis testing
parameters1 <- list()  # Parameters for pre-announcement period
parameters2 <- list()  # Parameters for post-announcement period
combined <- list()  # Combined results

# Parallel computing setup
totalCores <- detectCores() - 1  # Use all but one core
cluster <- makeCluster(totalCores) 
registerDoParallel(cluster)

# Main GPD parameter estimation and comparison process
Params_GPD <- foreach(i = 1:length(input_pre)) %dopar% {
  print(paste("Processing agent", i))
  
  x1n <- input_pre[[i]]  # Pre-announcement data for agent i
  x2n <- input_post[[i]]  # Post-announcement data for agent i
  
  if (length(x1n) > Z & length(x2n) > Z) {
    x1 <- scale(x1n, center = TRUE, scale = TRUE)  # Scaling pre-announcement data
    x2 <- scale(x2n, center = TRUE, scale = TRUE)  # Scaling post-announcement data
    
    # Compute Persistence Diagrams for pre and post announcement periods
    DiagGrid1[[i]] <- compute_persistence_diagram(x1, eta)
    DiagGrid2[[i]] <- compute_persistence_diagram(x2, eta)
    
    # GPD Parameter estimation and hypothesis testing for both periods
    params_result <- estimate_and_compare_gpd_parameters(DiagGrid1[[i]], DiagGrid2[[i]], K)
    
    # Store results
    pval[[i]] <- params_result$p_value
    parameters1[[i]] <- params_result$parameters1
    parameters2[[i]] <- params_result$parameters2
    combined[[i]] <- params_result$combined_matrix
  }
}

# Shutdown cluster after processing
stopCluster(cluster)

# ---------------------------------------------------------------------------------------------------
# Functions used in the main process
# ---------------------------------------------------------------------------------------------------

# Function to compute persistence diagram
compute_persistence_diagram <- function(data, eta) {
  Xlim <- range(data[, 1])
  Ylim <- range(data[, 2])
  by_val <- (Xlim[2] - Xlim[1]) / 1000
  Grid <- expand.grid(seq(from = Xlim[1], to = Xlim[2], by = by_val),
                      seq(from = Ylim[1], to = Ylim[2], by = by_val))
  diagram <- gridDiag(X = data, FUN = kde, h = eta, lim = cbind(Xlim, Ylim), by = by_val,
                      sublevel = FALSE, library = "Dionysus", printProgress = TRUE)
  return(diagram)
}

# Function to estimate GPD parameters and compare between two periods
estimate_and_compare_gpd_parameters <- function(DiagGrid1, DiagGrid2, K) {
  # Extract persistence data for GPD estimation
  X_updated1 <- preprocess_persistence_data(DiagGrid1$diagram)
  X_updated2 <- preprocess_persistence_data(DiagGrid2$diagram)
  
  # Estimate GPD parameters for pre-announcement period
  result1 <- estimate_gpd(X_updated1, K)
  optimal_parameters1 <- result1$parameters
  variances1 <- result1$variances
  
  # Estimate GPD parameters for post-announcement period
  result2 <- estimate_gpd(X_updated2, K)
  optimal_parameters2 <- result2$parameters
  variances2 <- result2$variances
  
  # Perform hypothesis testing between pre and post period parameters
  p_value <- perform_hypothesis_testing(optimal_parameters1, optimal_parameters2, variances1, variances2, K)
  
  # Combine results into a matrix
  combined_matrix <- rbind(p_value, optimal_parameters1, optimal_parameters2, variances1, variances2)
  
  return(list(p_value = p_value, parameters1 = optimal_parameters1, parameters2 = optimal_parameters2,
              combined_matrix = combined_matrix))
}

# Function to preprocess persistence data
preprocess_persistence_data <- function(diagram) {
  X_new <- diagram[which(diagram[, 1] == 0), c(2, 3)]
  X_updated <- X_new
  X_updated[, 1] <- X_new[, 2]
  X_updated[, 2] <- X_new[, 1]
  return(X_updated[-1, ])  # Remove the first entry and return
}

# Function to estimate GPD parameters
estimate_gpd <- function(data, K) {
  Hmat <- Hpi.diag(x = data)  # Bandwidth matrix for KDE
  fhat <- kde(x = data, H = as.matrix(Hmat))  # Kernel Density Estimation
  
  # Run GPD optimization
  optimal_func <- ucminf(par = initial_parameters, fn = neg_log_likelihood_mod, xn = data,
                         Hmat = Hmat, int_whole = compute_integral(data, Hmat), control = list(trace = 1, xtol = 1e-4, maxeval = 50), hessian = 3)
  
  optimal_parameters <- optimal_func$par
  optimal_hessinv <- optimal_func$invhessian
  
  return(list(parameters = optimal_parameters, variances = diag(optimal_hessinv)))
}

# Function to perform hypothesis testing between two parameter sets
perform_hypothesis_testing <- function(params1, params2, var1, var2, K) {
  test_statistic <- rep(0, K + 1)
  p_value <- rep(0, K + 1)
  
  for (j in 1:(K + 1)) {
    test_statistic[j] <- (params1[j] - params2[j]) / sqrt(var1[j] + var2[j])
    p_value[j] <- 2 * pnorm(abs(test_statistic[j]), lower.tail = FALSE)
  }
  
  return(p_value)
}
