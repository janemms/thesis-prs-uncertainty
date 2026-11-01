# ====
# Comparison of classifications given by Methods P and U
# ====

# Output: Saves classification metrics for both methods across all simulation 
#         replicates and parameter combinations.

library(here)
source(here("R", "thesis-functions.R"))
library(dplyr)

# ----
# 1. Simulation settings: define parameters and cost scenarios
# ----

n.target = 1000
nreps = 100

M.vals <- c(1000, 10000, 50000) #
N.vals = c(50000, 100000, 200000) # GWAS sample sizes
base.grid <- expand.grid(
  h2 = c(0.2, 0.5, 0.8), # heritabilities
  K = c(0.1, 0.2) # prevalences
)
cost.scenarios <- data.frame(
  C_FP = c(1, 1, 1, 1, 3, 5, 10),
  C_FN = c(10, 5, 3, 1, 1, 1, 1)
)

total.cells <- nrow(base.grid) * length(M.vals) * length(N.vals) * nrow(cost.scenarios)

cat("Total parameter cells:", total.cells, "\n")
cat("Replicates per cell:", nreps, "\n")
cat("Total iterations:", total.cells * nreps, "\n\n")

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
    true.effects <- simulate.true.effects(M = M, h2 = h2) # scaled causal effects
    genotypes <- simulate.genotypes(n = n.target, ps = true.effects$ps) 
    e <- rnorm(n.target, mean = 0, sd = sqrt(1 - h2))
    
    # Filter invariant SNPs and compute true.gv
    filtered.fixed <- filter.and.scale.geno(genotypes, true.effects$ps, true.betas = true.effects$true.betas) 
    X <- filtered.fixed$X
    true.betas <- filtered.fixed$true.betas
    idx.filter <- filtered.fixed$idx.filter
    true.gv <- as.numeric(true.betas %*% t(X))
    
    cat(sprintf("  -> SNPs kept: %d / %d\n", sum(idx.filter), length(idx.filter)))
    
    # liability and disease status
    g <- true.gv # genetic component
    liability <- g + e 
    t <- qnorm(1-K) # threshold for disease
    case.status <- as.integer(liability > t)
    
    for (n.idx in seq_along(N.vals)){
      N <- N.vals[n.idx]
      
      # ----
      # 3. Simulate GWAS data: gwas effect estimates, calculate posterior means and variances
      # ----
      for (rep in 1:nreps) {
        
        set.seed(get.seed(pop.row*1000 + n.idx*100 + M.idx * 10, rep))  # unique seed per row and rep
        
        # single GWAS draw: beta_hat ~ true_effects + sampling noise given N
        gwas <- simulate.gwas.stats(N, true.effects$ps, true.effects$true.betas, true.effects$var, K)
        
        # posterior mean and variance per snp
        post.mean <- (gwas$var.scaled / (gwas$var.scaled + gwas$sigma.est)) * gwas$beta.est
        post.var <- (gwas$var.scaled * gwas$sigma.est / (gwas$var.scaled + gwas$sigma.est)) 
        
        # filtered posterior mean and variance
        post.mean.filtered <- post.mean[idx.filter]
        post.var.filtered <- post.var[idx.filter]
        
        # compute prs with closed form posterior
        prs <- compute.prs.closed(X, post.mean.filtered, post.var.filtered, true.betas)
        
        # ----
        # 4. Predict disease risk and classify individuals
        # ----
        for (cost.row in 1:nrow(cost.scenarios)){
          
          C_FP <- cost.scenarios$C_FP[cost.row]
          C_FN <- cost.scenarios$C_FN[cost.row]
          
          prob.cutoff <- C_FP / (C_FP + C_FN)
          
          prob.over.P <- pnorm(t, mean = prs$pgs, sd = sqrt(1 - h2), lower.tail = FALSE) 
          prob.over.U <- pnorm(t, mean = prs$pgs, sd = sqrt(prs$pgs.var + 1 - h2), lower.tail = FALSE)
          
          hrisk.P <- which(prob.over.P > prob.cutoff)
          hrisk.U <- which(prob.over.U > prob.cutoff)
          
          # ----
          # 5. Compute metrics and save results
          # ----
          m.P <- compute.metrics(hrisk.P, case.status, C_FP, C_FN, n.target)
          m.U <- compute.metrics(hrisk.U, case.status, C_FP, C_FN, n.target)
          
          for (met in c("point.estimate", "uncertainty.aware") ) {
            
            if (met == "point.estimate"){
              m <- m.P
            } else {
              m <- m.U
            } 
            var(true.gv)
            var(liability)
            mean(case.status)
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

res.probabilistic.100rep.1 <- do.call(rbind, results)
saveRDS(res.probabilistic.100rep.1, here("results", "classification_sim_1.rds"))
