library(here)
source(here("R", "plot-functions.R"))

library(dplyr)
library(ggplot2)

# disease risk comparison between models M1 and M2

risks <- readRDS(here("results", "risk_comparison_all.rds"))
summary <- readRDS(here("results", "risk_comparison_summary.rds"))

# effect of N: fixed population, new gwas draw for each N

p.N <- risks %>% 
  filter(pop == 1, M == 10000, h2 == 0.8, K == 0.2) %>%
  mutate(panel = paste0("N = ", format(N, big.mark = ",", scientific = FALSE, trim = TRUE)),
         panel = factor(panel, levels = unique(panel))) %>% 
  plot.risk.scatter()

print(p.N)
ggsave(here("figures", "risk-comparison-N.png"),  p.N,  width = 10, height = 5, dpi = 300)

#risk.per.N %>% group_by(N) %>% summarise(mean.abs.diff = mean(abs(p.m2 - p.m1)), max.abs.diff = max(abs(p.m2 - p.m1)))


# effect of h2: independent populations for each h2
p.h2 <- risks %>% 
  filter(pop == 1, M == 10000, N == 50000, K == 0.2) %>%
  mutate(panel = paste0("h\u00B2 = ", h2),
         panel = factor(panel, levels = unique(panel))) %>% 
  plot.risk.scatter()

print(p.h2)
ggsave(here("figures", "risk-comparison-h2.png"), p.h2, width = 10, height = 5, dpi = 300)

# mean absolute differences 
#risk.per.h2 %>% group_by(h2) %>% summarise(mean.abs.diff = mean(abs(p.m2 - p.m1)), max.abs.diff = max(abs(p.m2 - p.m1)))


# effect of K: independent populations for each K
p.h2 <- risks %>% 
  filter(pop == 1, M == 10000, N == 50000, h2 == 0.5) %>%
  mutate(panel = paste0("K = ", format(K, big.mark = ",", scientific = FALSE, trim = TRUE)),
         panel = factor(panel, levels = unique(panel))) %>% 
  plot.risk.scatter()

print(p.h2)
