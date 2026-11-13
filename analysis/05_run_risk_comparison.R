
library(here)
source(here("R", "thesis-functions.R"))

# ==============================================================================
# Comparison of disease risks given by Methods P and U
#   Method P (point.estimate): risk from PRS only
#   Method U (uncertainty.aware): risk also includes posterior variance of genetic value
# Both methods use the same population and the same GWAS draw
# Risks are predicted under the liability threshold model
#
# Output: results/risk_comparison_all.rds   disease risks for both methods (one row per individual)
#         results/risk_comparison_summary.rds   mean and max absolute differences between risks give by methods
# ==============================================================================

# ----
# 1. Simulation settings: define parameters and cost scenarios
# ----
n.target <- 1000
n.pop <- 10
pop.grid <- expand.grid(M = c(1000, 10000), h2 = c(0.2, 0.5, 0.8))
K.vals <- c(0.1, 0.2)
N.vals <- c(50000, 200000)

risk.list <- list()
summary.list <- list()

for (grid.row in 1:nrow(pop.grid)) {
  M <- pop.grid$M[grid.row]
  h2 <- pop.grid$h2[grid.row]
  
  for (pop.idx in 1:n.pop) {
    set.seed(get.seed(grid.row, pop.idx))
    # ----
    # 2. Simulate the target population: causal effects, genotypes, scaling
    # ----
    pop <- simulate.population(M, h2, n.target)
    
    for (K in K.vals) {
      for (N in N.vals) {
        # ----
        # 3. Simulate GWAS, compute PRS and predict disease risks under LTM for both methods
        # ----
        prs <- prs.posterior.var(pop, N, K)
        probs <- disease.probs(prs, h2, K)
        difference <- abs(probs$p.mU - probs$p.mP)
        
        # ----
        # 4. Compute metrics and save results
        # ----
        summary.list[[length(summary.list) + 1]] <- data.frame(M, h2, K, N, pop.idx = pop.idx, mean.abs.diff = mean(difference), max.abs.diff = max(difference))
        risk.list[[length(risk.list) + 1]] <- cbind(M = M, h2 = h2, K = K, N = N, pop.idx = pop.idx, id = 1:n.target, round(probs, 4))
      }
      
    }
  }
}

# save disease probabilities and summary statistics to rds files
saveRDS(do.call(rbind, risk.list), here("results", "risk_comparison_all.rds"), compress = "xz")
saveRDS(do.call(rbind, summary.list), here("results", "risk_comparison_summary.rds"))
