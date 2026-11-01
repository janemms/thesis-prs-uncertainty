# comparison of disease risks given by two PRS classification methods 
library(here)
source(here("R", "thesis-functions.R"))


n.target <- 1000
n.pop <- 10
pop.grid <- expand.grid(M = c(1000, 10000), h2 = c(0.2, 0.5, 0.8))
K.vals <- c(0.1, 0.2)
N.vals <- c(50000, 200000)

risk.list <- list()
summary.list <- list()

for (i in 1:nrow(pop.grid)) {
  M <- pop.grid$M[i]
  h2 <- pop.grid$h2[i]
  
  for (p in 1:n.pop) {
    set.seed(get.seed(i, p))
    pop <- simulate.population(M, h2, n.target)
    
    for (K in K.vals) {
      for (N in N.vals) {
        prs <- prs.posterior.var(pop, N, K)
        probs <- disease.probs(prs, h2, K)
        difference <- abs(probs$p.m2 - probs$p.m1)
        
        summary.list[[length(summary.list) + 1]] <- data.frame(M, h2, K, N, pop = p, mean.abs.diff = mean(difference), max.abs.diff = max(difference))
        risk.list[[length(risk.list) + 1]] <- cbind(M = M, h2 = h2, K = K, N = N, pop = p, id = 1:n.target, round(probs, 4))
      }
      
    }
  }
}

saveRDS(do.call(rbind, risk.list), here("results", "risk_comparison_all.rds"), compress = "xz")
saveRDS(do.call(rbind, summary.list), here("results", "risk_comparison_summary.rds"))
