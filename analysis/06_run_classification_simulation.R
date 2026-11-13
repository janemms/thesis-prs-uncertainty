
library(here)
source(here("R", "thesis-functions.R"))
library(dplyr)

# ==============================================================================
# Comparison of classifications given by Methods P and U
#   Method P (point.estimate): risk from PRS only
#   Method U (uncertainty.aware): risk also includes posterior variance of genetic value
# Both methods use the same population and the same GWAS draw
#
# Output: results/classification_sim_seed.rds
#         Classification metrics for both methods across all simulation replicates
#         and parameter combinations.
# ==============================================================================

# ----
# 1. Simulation settings: define parameters and cost scenarios
# ----

n.target <- 1000
nreps <- 100

M.vals <- c(1000, 10000, 50000) # number of SNPs
N.vals = c(50000, 100000, 200000) # GWAS sample sizes
base.grid <- expand.grid(
  h2 = c(0.2, 0.5, 0.8), # heritabilities
  K = c(0.1, 0.2) # prevalences
)
cost.scenarios <- data.frame(
  C_FP = c(1, 1, 1, 1, 3, 5, 10), # cost of a false positive
  C_FN = c(10, 5, 3, 1, 1, 1, 1) # cost of a false negative
)

risk.cols <- c(point.estimate = "p.mP", uncertainty.aware = "p.mU")

total.cells <- nrow(base.grid) * length(M.vals) * length(N.vals) * nrow(cost.scenarios)

cat("Total parameter cells:", total.cells, "\n")
cat("Replicates per cell:", nreps, "\n")
cat("Total iterations:", total.cells * nreps, "\n\n")

# one entry per iteration and method
results <- vector("list", total.cells * nreps * 2)
res.idx <- 0

start.time <- Sys.time()
for (M.idx in seq_along(M.vals)) {
  M <- M.vals[M.idx]
  
  for (pop.row in 1:nrow(base.grid)) {
    h2 <- base.grid$h2[pop.row]
    K <- base.grid$K[pop.row]
    
    set.seed(get.seed(pop.row * 10 + M.idx, 1)) # unique seed per fixed population
    
    elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "mins"))
    cat(sprintf("[Cell %d/%d] h2=%.2f K=%.2f | elapsed: %.1f min",
                pop.row, nrow(base.grid), h2, K, elapsed))
    
    # ----
    # 2. Simulate the target population: causal effects, genotypes, disease statuses
    # ----
    pop <- simulate.population(M, h2, n.target)
    case.status <- simulate.case.status(pop, h2, K)
    
    for (n.idx in seq_along(N.vals)){
      N <- N.vals[n.idx]
      
      # ----
      # 3. Simulate GWAS data: draw GWAS effect estimates, calculate posterior means and variances
      # ----
      for (rep in 1:nreps) {
        
        set.seed(get.seed(pop.row*1000 + n.idx*100 + M.idx * 10, rep))  # unique seed per row and rep
        
        prs <- prs.posterior.var(pop, N, K)
                
        # ----
        # 4. Predict disease risk and classify individuals
        #     Method P: liability | PRS ~ N(PRS, 1 - h2)
        #     Method P: liability | gwas data ~ N(PRS, 1 - h2 + prs.var)
        # ----
        
        probs <- disease.probs(prs, h2, K)
        
        for (cost.row in 1:nrow(cost.scenarios)){
          
          C_FP <- cost.scenarios$C_FP[cost.row]
          C_FN <- cost.scenarios$C_FN[cost.row]
          
          prob.cutoff <- C_FP / (C_FP + C_FN)
      
          # ----
          # 5. Classify, compute metrics and save results
          # ----
       
          for (met in c("point.estimate", "uncertainty.aware") ) {
            
            hrisk <- which(probs[[risk.cols[[met]]]] > prob.cutoff)
            m <- compute.metrics(hrisk, case.status, C_FP, C_FN, n.target)

            res.idx <- res.idx + 1
            results[[res.idx]] <- data.frame(
              rep = rep, N = N, M = M, n = n.target, h2 = h2, K = K,
              cost.FP = C_FP, cost.FN = C_FN,
              method = met,
              m
            )
          }
        }
      }
    }
  }
}
cat("\nTotal time:", as.numeric(difftime(Sys.time(), start.time, units = "mins")), "minutes\n")

res.classification.1 <- do.call(rbind, results)
# save to rds file, name: "classification_sim_seed"
saveRDS(res.classification.1, here("results", "classification_sim_1.rds"))
