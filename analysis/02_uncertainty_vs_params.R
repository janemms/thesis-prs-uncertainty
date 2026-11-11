library(here)
source(here("R", "thesis-functions.R"))
source(here("R", "plot-functions.R"))

library(ggplot2)
library(dplyr)
library(tidyr)

# ==============================================================================
# Effect of N, M and h2 on the posterior variance of genetic value
# Output: figures/var-vs-N.png, figures/var-vs-M.png, figures/var-vs-h2.png
# ==============================================================================

# fixed settings used on all plots
n.target <- 1000
K <- 0.2

# default parameter values
M.def <- 10000
h2.def <- 0.5
N.def <- 100000

# format numbers with thousand separators for axis labels
add.sep <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")

# ----
# 1. Effect of GWAS sample size N 
#     Fixed: one target population (M, h2), Varied: N, new GWAS draw for each N
# ----
set.seed(111)
N.vals <- c(50000, 100000, 200000, 300000, 400000)
pop <- simulate.population(M = M.def, h2 = h2.def, n.target)
var.vs.N <- sweep.var(values = N.vals, function(N) prs.posterior.var(pop, N, K)$prs.var)
plot.N <- plot.violin(var.vs.N, "N", add.sep) 
print(plot.N)

ggsave(here("figures", "var-vs-N.png"), plot = plot.N, width = 8, height = 8, dpi = 300)

# ----
# 2. Effect of number of SNPs M 
#     Fixed: h2, N, Varied: M, new population (and GWAS draw) for each M
# ----
set.seed(111)
M.vals <- c(1000, 10000, 50000, 100000)
var.vs.M <- sweep.var(M.vals, function(M) prs.posterior.var(simulate.population(M, h2 = h2.def, n.target), N.def, K)$prs.var)
plot.M <- plot.violin(var.vs.M, "M", add.sep) 
print(plot.M)

ggsave(here("figures", "var-vs-M.png"), plot = plot.M, width = 8, height = 8, dpi = 300)

# ----
# 3. Effect of heritability h2
#     Fixed: M, N, Varied: h2, new population (and GWAS draw) for each h2
# ----

set.seed(111)
h2.vals <- c(0.1, 0.2, 0.4, 0.6, 0.8)
var.vs.h2 <- sweep.var(h2.vals, function(h2) prs.posterior.var(simulate.population(M.def, h2, n.target), N.def, K)$prs.var)
plot.h2 <- plot.violin(var.vs.h2, expression(h^2))
print(plot.h2)

ggsave(here("figures", "var-vs-h2.png"), plot = plot.h2, width = 8, height = 8, dpi = 300)
