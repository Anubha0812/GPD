
# GPD Parameter Estimation and Hypothesis Testing
This code can be used to compare the distributions of any two persistence diagrams.
This repository contains R code to compare compares the Gibbs Distribution parameter estimates of the persistence diagrams for transaction data between pre-announcement and post-announcement periods for multiple agents based on their transaction data. The code performs the following tasks:

- Estimation of GPD parameters for each agent in both periods.
- Hypothesis testing between parameters from the pre-announcement and post-announcement periods.
- Output of p-values, parameter estimates, and variances in a combined matrix.

## Table of Contents

- [Requirements](#requirements)
- [Input Data](#input-data)
- [Code Structure](#code-structure)
- [Usage](#usage)
- [Output](#output)

---

## Requirements

To run the code, you'll need to have the following installed on your system:

- R (version 3.6+)
- Required R packages:
  - `foreach`
  - `doParallel`
  - `parallel`
  - `TDA`
  - `pracma`
  - `dbscan`
  - `ucminf`
  - `mvtnorm`
  - `ks`

You can install the required packages using the following commands in R:

```r
install.packages(c("foreach", "doParallel", "parallel", "TDA", "pracma", "dbscan", "ucminf", "mvtnorm", "ks"))
```

---

## Input Data

The input data consists of two lists of transaction data:

1. `input_pre`: A list where each element contains transaction data for an agent in the pre-announcement period.
2. `input_post`: A list where each element contains transaction data for an agent in the post-announcement period.

### Format:

- Each element of `input_pre` and `input_post` should be a matrix or data frame where rows represent transactions, and columns represent relevant features of the transaction.

### Example:

```r
input_pre <- list(agent1_data, agent2_data, ...)
input_post <- list(agent1_data, agent2_data, ...)
```
---

## Code Structure

- **Main Script**: The script iterates over all agents in the input data, estimates GPD parameters for both the pre-announcement and post-announcement periods, and conducts hypothesis testing.
  
- **Functions**: 
  - `compute_persistence_diagram`: Computes persistence diagrams for the transaction data.
  - `estimate_gpd`: Estimates the GPD parameters for each period using a specified optimization procedure.
  - `perform_hypothesis_testing`: Performs statistical tests to compare the GPD parameters between pre and post-announcement periods.

- **Parallelization**: The code uses parallel processing to distribute the computation across multiple cores for faster execution.

---


## Output

For each agent, the following results are generated:

1. **Parameter Estimates**: GPD parameter estimates for the pre and post-announcement periods.
2. **P-values**: Hypothesis test p-values for comparing parameter estimates between periods.
3. **Combined Matrix**: A matrix for each agent containing:
   - p-values,
   - GPD parameter estimates (for both periods),
   - Variances of the estimates.

You can access the results in the variables `parameters1`, `parameters2`, `pval`, and `combined`.

---

