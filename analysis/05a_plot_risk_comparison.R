library(here)
source(here("R", "plot-functions.R"))

library(dplyr)
library(ggplot2)

# ==============================================================================
# Comparison of disease risks predicted by Method P (point estimate) and 
#                                          Method U (uncertainty-aware).
# Input: results/risk_comparison_all.rds
# Output: figures/risk-comparison-N.png, figures/risk-comparison-h2.png
# ==============================================================================

# simularion results: one row per individual and parameter setting
risks <- readRDS(here("results", "risk_comparison_all.rds"))

# ----
# 1. Effect of GWAS sample size N
#     Fixed: population, M, h2, K, Varied: N, new GWAS draw for each N
# ----
p.N <- risks %>% 
  filter(pop.idx == 1, M == 10000, h2 == 0.8, K == 0.2) %>%
  mutate(panel = paste0("N = ", format(N, big.mark = ",", scientific = FALSE, trim = TRUE)),
         panel = factor(panel, levels = unique(panel))) %>% 
  plot.risk.scatter()

print(p.N)
ggsave(here("figures", "risk-comparison-N.png"),  p.N,  width = 10, height = 5, dpi = 300)

# ----
# 2. Effect of heritability h2
#     Fixed: M, N, K, Varied: h2, new population (and GWAS draw) for each h2
# ----
p.h2 <- risks %>% 
  filter(pop.idx == 1, M == 10000, N == 50000, K == 0.2) %>%
  mutate(panel = paste0("h\u00B2 = ", h2),
         panel = factor(panel, levels = unique(panel))) %>% 
  plot.risk.scatter()

print(p.h2)
ggsave(here("figures", "risk-comparison-h2.png"), p.h2, width = 10, height = 5, dpi = 300)

