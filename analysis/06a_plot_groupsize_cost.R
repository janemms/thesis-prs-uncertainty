library(here)
source(here("R", "thesis-functions.R"))
source(here("R", "plot-functions.R"))

# ==============================================================================
# Comparison of classifications given by Methods P and U across cost scenarios
#   Boxplots over simulation replicates of
#     - the proportion of individuals classified as high risk
#     - the realized cost per individual
#   for different GWAS sample sizes N and heritabilities h2.
#
# Input: results/classification_sim_seed.rds
# Output: figures/prop-selected-N.png, figures/prop-selected-h2.png,
#         figures/realized-cost-N.png, figures/realized-cost-h2.png
# ==============================================================================

# one row per replicate, parameter combination, cost scenario and method
results.df <- readRDS(here("results", "classification_sim_1.rds"))

# ----
# 1. Proportion classified as high risk
# ----

# effect of N
prop.N <- results.df %>%
  prepare.results(M = 10000, K = 0.2, N = c(50000, 200000), h = 0.8) %>% 
  plot.comparison.boxplot("N", "prop.selected", "Proportion classified as high risk")
print(prop.N)
ggsave(here("figures", "prop-selected-N.png"), plot = prop.N, width = 8, height = 5, dpi = 300)

# effect of h2
prop.h2 <- results.df %>%
  prepare.results(M = 10000, K = 0.2, N = 50000, h = c(0.2, 0.8)) %>% 
  plot.comparison.boxplot("h2", "prop.selected", "Proportion classified as high risk")
print(prop.h2)
ggsave(here("figures", "prop-selected-h2.png"), plot = prop.h2, width = 8, height = 5, dpi = 300)

# ----
# 2. Realized cost per individual
# ----

# effect of N
cost.N <- results.df %>%
  prepare.results(M = 10000, K = 0.2, N = c(50000, 200000), h = 0.5) %>% 
  plot.comparison.boxplot("N", "cost.per.ind", "Realized cost")
print(cost.N)
ggsave(here("figures", "realized-cost-N.png"), plot = cost.N, width = 8, height = 5, dpi = 300)

# effect of h2
cost.h2 <- results.df %>%
  prepare.results(M = 10000, K = 0.2, N = 50000, h = c(0.2, 0.8)) %>% 
  plot.comparison.boxplot("h2", "cost.per.ind", "Realized cost")
print(cost.h2)
ggsave(here("figures", "realized-cost-h2.png"), plot = cost.h2, width = 8, height = 5, dpi = 300)
