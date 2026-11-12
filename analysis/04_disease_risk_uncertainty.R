library(here)
library(ggplot2)
library(dplyr)
source(here("R", "thesis-functions.R"))

# ==============================================================================
# Uncertainty in PRS propagated to disease risk
# Output: figures/disease-risk-uncertainty.png
# ==============================================================================

# ----
# 1. Settings
# ----
M = 1000
n.target = 1000
N = 200000
h2 = 0.5
Ks = c(0.1, 0.2)

# ----
# 2. Disease risk with 95% credible interval, function of the PRS
# ----
# risk at the posterior mean genetic value and at the bounds of its credible interval
#   prs: output of prs.posterior.var (PRS and posterior variance of genetic value)
#   K: disease prevalence
#   h2: heritability
#   level: credible level, default 95%
# Returns:
#   data frame with prs, risk, risk.low, risk.up and K
risk.interval <- function(prs, K, h2, level = 0.95) {
  t <- qnorm(1 - K) # liability threshold
  z <- qnorm(1 - (1 - level)/2)
  # P(liability > t | gv), residual variance 1 - h2
  risk <- function(gv) pnorm( (gv - t) / (sqrt(1-h2)) )
  data.frame(prs = prs$prs, 
             risk = risk(prs$prs),
             risk.low = risk(prs$prs - z * sqrt(prs$prs.var)),
             risk.up = risk(prs$prs + z * sqrt(prs$prs.var)),
             K = factor(K))
}

# one population, a new GWAS draw for each K
set.seed(1)
pop <- simulate.population(M, h2, n.target)

plot.df <- do.call(rbind, lapply(Ks, function(K){
  risk.interval(prs.posterior.var(pop, N, K), K, h2)
}))

# ----
# 3. Plot: risk curve with creadible inteval band, one color per prevalence
# ----
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
  theme(panel.border = element_blank(), 
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5), 
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(plot.risk.interval)


ggsave(here("figures", "disease-risk-uncertainty.png"), plot = plot.risk.interval, width = 8, height = 5, dpi = 300)
