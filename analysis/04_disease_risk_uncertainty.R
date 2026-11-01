library(here)
source(here("R", "thesis-functions.R"))


# ------------------------------------------------------------------------------
# ----------------- uncertainty propagated to disease risk ---------------------
# ------------------------------------------------------------------------------
M = 1000
n.target = 1000
N = 200000
h2 = 0.5
Ks = c(0.1, 0.2)


# risk at posterior mean GV and at the bounds of its 95 % credible interval
risk.interval <- function(prs, K, h2, level = 0.95) {
  t <- qnorm(1 - K)
  z <- qnorm(1 - (1 - level)/2)
  risk <- function(gv) pnorm( (gv - t) / (sqrt(1-h2)) )
  data.frame(prs = prs$pgs, 
             risk = risk(prs$pgs),
             risk.low = risk(prs$pgs - z * sqrt(prs$pgs.var)),
             risk.up = risk(prs$pgs + z * sqrt(prs$pgs.var)),
             K = factor(K))
}

set.seed(1)
pop <- simulate.population(M, h2, n.target)

plot.df <- do.call(rbind, lapply(Ks, function(K){
  risk.interval(prs.posterior.var(pop, N, K), K, h2)
}))

ribbon.col <- c("0.1" = "#CC79A7", "0.2" = "#56B4E9" )
line.col <- c("0.1" = "red4", "0.2" = "darkblue" )

plot.risk.interval <- ggplot(plot.df, aes(x = prs)) +
  geom_ribbon(aes(ymin = risk.low, ymax = risk.up, fill = K), alpha = 0.4, color = NA) +
  geom_line(aes(y = risk, color = K), linewidth = 1) +
  scale_color_manual(values = line.col) +
  scale_fill_manual(values = ribbon.col) +
  labs(x = "PRS",
    y = "Estimated disease risk",
    color = "Prevalence K",
    fill = "Prevalence K"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(plot.risk.interval)


ggsave(here("figures", "disease-risk-uncertainty.png"), plot = plot.risk.interval, width = 8, height = 5, dpi = 300)



true.effects <- simulate.true.effects(M = M, h2 = h2) # scaled causal effects
genotypes <- simulate.genotypes(n = n.target, ps = true.effects$ps)
e <- rnorm(n.target, mean = 0, sd = sqrt(1 - h2))

# Filter invariant SNPs and compute true.gv
filtered.fixed <- filter.and.scale.geno(genotypes, ps = true.effects$ps, true.betas = true.effects$true.betas) 
X <- filtered.fixed$X
true.betas <- filtered.fixed$true.betas
idx.filter <- filtered.fixed$idx.filter

true.gv <- as.numeric(true.betas %*% t(X))
liability <- true.gv + e 

plot.list <- vector("list", length(Ks))

for (i in seq_along(Ks)) {
  K <- Ks[i]

  t <- qnorm(1-K) # threshold for disease
  case.status <- as.integer(liability > t)
  
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
  
  q.low <- 0.025
  q.up <- 0.975
  prs.low <- prs$pgs + sqrt(prs$pgs.var) * qnorm(q.low)
  prs.up <- prs$pgs + sqrt(prs$pgs.var) * qnorm(q.up)
  
  risk.low <- pnorm((prs.low - t)/sqrt(1-h2))
  risk.up <- pnorm((prs.up - t)/sqrt(1-h2))
  risk.median <- pnorm((prs$pgs - t)/sqrt(1-h2))
  
  plot.list[[i]] <- data.frame(
    prs = prs$pgs,
    risk.median = risk.median,
    risk.low = risk.low,
    risk.up = risk.up,
    case.status = factor(case.status, labels = c("Unaffected", "Affected")), 
    K = paste0(K)
  ) %>%
    arrange(prs) %>%
    mutate(prs.quantile = row_number() / n())  
}

plot.df.all <- do.call(rbind, plot.list)
plot.df.all$K <- factor(plot.df.all$K)

ribbon.col <- c("0.1" = "CC79A7#", "0.2" = "#56B4E9" )
line.col <- c("0.1" = "red4", "0.2" = "darkblue" )


