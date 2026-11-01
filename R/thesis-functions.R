# ----------------------------------------------------------------------------------------------------
# Function for drawing genetic architecture (allele frequencies + true causal effect sizes) 
# for a given heritability h2 and number of SNPs M. 
# ----------------------------------------------------------------------------------------------------
simulate.true.effects <- function(M, h2, ps = NULL){
  if (is.null(ps)) {
    ps <- runif(M, 0.01, 0.99)
  }
  
  var <- h2 / ((2 * ps * (1 - ps)) * M)
  
  # draw true betas with variance var, then scale
  betas <- rnorm(M, mean = 0, sd = sqrt(var))
  betas <- betas * sqrt(2 * ps * (1 - ps))
  
  return(list(
    ps = ps, # allele frequencies
    true.betas = betas, # true (scaled) causal effect sizes
    var = var # unscaled per-SNP variance
  ))
}

# ----------------------------------------------------------------------------------------------------
# Function for generating genotypes for n individuals given allele frequencies ps.
# ----------------------------------------------------------------------------------------------------
simulate.genotypes <- function(n, ps){
  M <- length(ps)
  X <- matrix(rbinom(n * M, size = 2, prob = rep(ps, each = n)), nrow = n, ncol = M)
  return(X)
}

# ----------------------------------------------------------------------------------------------------
# Function for generating GWAS summary statistics for a fixed genetic architecture.
# ----------------------------------------------------------------------------------------------------

simulate.gwas.stats <- function(N, ps, true.betas, var, K){
  # N: GWAS sample size
  # ps: allele frequencies (from simulate.true.effects)
  # true.betas: true causal effect sizes (from simulate.true.effects)
  # var: unscaled per-SNP variance (from simulate.true.effects)
  
  t <- qnorm(1-K)
  M <- length(ps)
  var.scaled <- var * 2 * ps * (1 - ps)
  var.scaled <- pmin(var.scaled, 0.999999)
  
  scaling.fact <- K*(1-K)/((dnorm(t))^2) 
  
  sigma.est <- scaling.fact/N
  sigma.est <- pmax(sigma.est, 1e-8)
  
  # single GWAS draw: beta_hat ~ N(true_beta, sigma.est)
  beta.est <- rnorm(M, mean = true.betas, sd = sqrt(sigma.est))
  
  return(list(
    beta.est = beta.est, 
    sigma.est = sigma.est, 
    var.scaled = var.scaled,
    ps = ps
  ))
}


# ----------------------------------------------------------------------------------------------------
# Function for calculating posterior mean (PRS) and variance of genetic value
# ----------------------------------------------------------------------------------------------------
compute.prs.closed <- function(X, post.mean, post.var, true.betas) {
  pgs.mean <- as.numeric(post.mean %*% t(X))
  pgs.var  <- as.numeric((X^2) %*% post.var)
  #true.gv  <- as.numeric(true.betas %*% t(X))
  list(pgs = pgs.mean, pgs.var = pgs.var)#, true.gv = true.gv)
}


# ----------------------------------------------------------------------------------------------------
# Function for calculating summary metrics for a classification given classification costs
# ----------------------------------------------------------------------------------------------------
compute.metrics <- function(hrisk, case.status, C_FP, C_FN, n) {
  selected <- logical(n)
  selected[hrisk] <- TRUE
  
  FP  <- sum( selected & case.status == 0)
  FN  <- sum(!selected &  case.status == 1)
  TP  <- sum( selected &  case.status == 1)
  TN  <- sum(!selected & case.status == 0)
  PPV <- if (sum(selected) > 0) TP / sum(selected) else NA
  NPV <- if (sum(!selected) > 0) TN / sum(!selected) else NA
  TPR <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  
  cost <- FP * C_FP + FN * C_FN   # realized cost
  
  data.frame(FP = FP, FN = FN, TP = TP, TN = TN,
             PPV = PPV, NPV = NPV, TPR = TPR,
             cost = cost,   
             n.selected = sum(selected))
}

# ----------------------------------------------------------------------------------------------------
# Function for filtering out invariant SNPs and standardizing genotypes.
# Returns filtered X plus correspondingly filtered architecture/GWAS/posterior objects.
# ----------------------------------------------------------------------------------------------------

filter.and.scale.geno <- function(X, ps, ...){
  
  idx.filter <- apply(X, 2, sd) > 0 # filter our SNPs with no variation
  X <- X[, idx.filter]
  ps.filt <- ps[idx.filter] 
  X <- sweep(X, 2, 2 * ps.filt, "-")
  X <- sweep(X, 2, sqrt(2 * ps.filt * (1-ps.filt)), "/")
  
  extras <- list(...)
  extras.filtered <- lapply(extras, function(v) {
    if (is.matrix(v)){
      v[, idx.filter, drop = FALSE]
    } else {
      v[idx.filter]
    }
  })
  
  return(c(list(X = X, idx.filter = idx.filter), extras.filtered))
}

# simulate one target population: effects, genotypes, filtering, scaling
simulate.population <- function(M, h2, n.target) {
  effects <- simulate.true.effects(M = M, h2 = h2)
  geno <- simulate.genotypes(n.target, ps = effects$ps)
  filt <- filter.and.scale.geno(geno, effects$ps, true.betas = effects$true.betas)
  
  list(effects = effects, X = filt$X, true.betas = filt$true.betas, idx = filt$idx.filter)
}

# simulate disease status for a population
simulate.case.status <- function(pop, h2, K) {
  true.gv <- as.numeric(pop$true.betas %*% t(pop$X))
  e <- rnorm(length(true.gv), mean = 0, sd = sqrt(1 - h2))
  liability <- true.gv + e
  as.integer(liability > qnorm(1-K))
}

# compute posterior variance of gv for each individual for a single GWAS draw
prs.posterior.var <- function(pop, N, K) {
  e <- pop$effects
  gwas <- simulate.gwas.stats(N, e$ps, e$true.betas, e$var, K)
  post.mean <- (gwas$var.scaled / (gwas$var.scaled + gwas$sigma.est)) * gwas$beta.est
  post.var <- (gwas$var.scaled * gwas$sigma.est / (gwas$var.scaled + gwas$sigma.est)) 
  compute.prs.closed(pop$X, post.mean[pop$idx], post.var[pop$idx], pop$true.betas)
}

# compute disease risks under the two methods
disease.probs <- function(prs, h2, K) {
  t <- qnorm(1 - K)
  data.frame(
    p.m1 = pnorm(t, mean = prs$pgs, sd = sqrt(1 - h2), lower.tail = FALSE),
    p.m2 = pnorm(t, mean = prs$pgs, sd = sqrt(1 - h2 + prs$pgs.var), lower.tail = FALSE)
  )
}

# expected posterior variance of GV
# expected.post.var <- function(M, h2, N, K) {
#   num <- K*(1-K)
#   phi.t.2 <- dnorm(qnorm(1-K))^2
#   denom <- num/h2 + N*phi.t.2/M
#   num/denom
# }


# get unique seed for a (row, rep) pair
get.seed <- function(row, rep, base.seed = 2025071) {
  base.seed + row * 100000L + rep
}
