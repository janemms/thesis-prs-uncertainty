# ------------------------------------------------------------------------------
# Draw genetic architecture (allele frequencies + true causal effect sizes). 
#   M: number of SNPs
#   h2: heritability
#   ps: allele frequencies (drawn if ps = NULL)
# Returns:
#   ps: allele frequencies
#   true.betas: scaled causal genetic effects
#   var: variance per SNP (unscaled)    
# ------------------------------------------------------------------------------
simulate.true.effects <- function(M, h2, ps = NULL){
  if (is.null(ps)) {
    ps <- runif(M, 0.01, 0.99)
  }
  
  var <- h2 / ((2 * ps * (1 - ps)) * M)
  
  # draw true betas with variance var, then scale
  betas <- rnorm(M, mean = 0, sd = sqrt(var))
  betas <- betas * sqrt(2 * ps * (1 - ps))
  
  return(list(
    ps = ps, 
    true.betas = betas, 
    var = var 
  ))
}

# ------------------------------------------------------------------------------
# Draw genotypes (0,1 or 2 effects alleles) for n individuals.
#   n: # of individuals
#   ps: allele frequencies
# Returns:
#   X: n x M genotype matrix 
# ------------------------------------------------------------------------------
simulate.genotypes <- function(n, ps){
  M <- length(ps)
  X <- matrix(rbinom(n * M, size = 2, prob = rep(ps, each = n)), nrow = n, ncol = M)
  return(X)
}

# ------------------------------------------------------------------------------
# Draw GWAS summary statistics for a fixed genetic architecture.
#   N: GWAS sample size
#   ps: allele frequencies 
#   true.betas: causal effect sizes 
#   var: variance per SNP (unscaled)    
#   K: disease prevalence
# Returns:
#   beta.est: estimated genetic effects from GWAS
#   sigma.est: sampling variance of the estimates
#   var.scaled: prior variance per SNP (scaled)
#   ps: allele frequencies
# ------------------------------------------------------------------------------
simulate.gwas.stats <- function(N, ps, true.betas, var, K){
   
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

# ------------------------------------------------------------------------------
# Calculate posterior mean (PRS) and variance of genetic value
#   X: n x M genotype matrix (standardized)
#   post.mean: posterior mean of each SNP effect
#   post.var: posterior variance of each SNP effect
# Returns:
#   prs: PRS of each individual
#   prs.var: posterior variance of genetic value for each individual
# ------------------------------------------------------------------------------
compute.prs <- function(X, post.mean, post.var) {
  prs.mean <- as.numeric(post.mean %*% t(X))
  prs.var  <- as.numeric((X^2) %*% post.var)
  list(prs = prs.mean, prs.var = prs.var)
}


# ------------------------------------------------------------------------------
# Calculate classification metrics and realized costs for a high-risk selection.
#   hrisk: indices of individuals classified as high risk
#   case.status: true 0/1 disease status of all individuals
#   C_FP: cost of a false positive classification
#   C_FN: cost of a false negative classification
#   n: # of classified individuals
# Returns:
#   data frame with FP, FN, TP, TN, PPV, NPV, TPR, cost, n.selected
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Remove invariant SNPs and standardize genotypes.
#   X: genotype matrix
#   ps: allele frequencies
#   ...: other vectors/matrices to be filtered
# Returns:
#   X: filtered and scaled genotype matrix
#   idx.filter: logical vector of kept SNPs
#   ...: other filtered objects
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# Simulate one target population: effects, genotypes, filtering, scaling.
#   M: # of SNPs
#   h2: heritability
#   n.target: # of individuals
# Returns:
#   effects: output of simulate.true.effects (ps, true.betas, var)
#   X: scaled genotypes
#   true.betas: filtered causal effects
#   idx: logical vector of kept SNPs
# ------------------------------------------------------------------------------
simulate.population <- function(M, h2, n.target) {
  effects <- simulate.true.effects(M = M, h2 = h2)
  geno <- simulate.genotypes(n.target, ps = effects$ps)
  filt <- filter.and.scale.geno(geno, effects$ps, true.betas = effects$true.betas)
  
  list(effects = effects, X = filt$X, true.betas = filt$true.betas, idx = filt$idx.filter)
}

# ------------------------------------------------------------------------------
# Simulate disease status from the liability threshold model.
#   pop: output of simulate.population
#   h2: heritability
#   K: disease prevalence
# Returns:
#   a vector of 0/1 disease statuses
# ------------------------------------------------------------------------------
simulate.case.status <- function(pop, h2, K) {
  true.gv <- as.numeric(pop$true.betas %*% t(pop$X))
  e <- rnorm(length(true.gv), mean = 0, sd = sqrt(1 - h2))
  liability <- true.gv + e
  as.integer(liability > qnorm(1-K))
}

# ------------------------------------------------------------------------------
# Calculate each individual's PRS and the posterior variance of genetic value 
# for a single GWAS draw.
#   pop: output of simulate.population
#   N: GWAS sample size
#   K: disease prevalence
# Returns:
#   prs: PRS of each individual
#   prs.var: posterior variance of genetic value per individual
# ------------------------------------------------------------------------------
prs.posterior.var <- function(pop, N, K) {
  e <- pop$effects
  gwas <- simulate.gwas.stats(N, e$ps, e$true.betas, e$var, K)
  post.mean <- (gwas$var.scaled / (gwas$var.scaled + gwas$sigma.est)) * gwas$beta.est
  post.var <- (gwas$var.scaled * gwas$sigma.est / (gwas$var.scaled + gwas$sigma.est)) 
  compute.prs(pop$X, post.mean[pop$idx], post.var[pop$idx])
}

# ------------------------------------------------------------------------------
# Compute disease risk for each individual under two methods.
#   prs: output of compute.prs (prs, prs.var)
#   h2: heritability
#   K: disease prevalence
# Returns:
#   p.mP: risk under LTM using PRS only
#   p.mU: risk under LTM when including posterior variance of genetic value
# ------------------------------------------------------------------------------
disease.probs <- function(prs, h2, K) {
  t <- qnorm(1 - K)
  data.frame(
    p.mP = pnorm(t, mean = prs$prs, sd = sqrt(1 - h2), lower.tail = FALSE),
    p.mU = pnorm(t, mean = prs$prs, sd = sqrt(1 - h2 + prs$prs.var), lower.tail = FALSE)
  )
}

# ------------------------------------------------------------------------------
# Draw a unique, reproducible seed for a (row, rep) pair.
#   row: row index of simulation grid
#   rep: replicate number
#   base.seed: initial seed of the simulation 
# Returns:
#   integer seed for R Random RNG
get.seed <- function(row, rep, base.seed = 2025071) {
  base.seed + row * 100000L + rep
}
