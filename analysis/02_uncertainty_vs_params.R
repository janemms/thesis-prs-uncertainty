library(here)
source(here("R", "thesis-functions.R"))
source(here("R", "plot-functions.R"))

library(ggplot2)
library(dplyr)
library(tidyr)

n.target = 1000
K <- 0.2
add.sep <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")

# Effect of N on prs var: one fixed population, new gwas for each N
set.seed(111)
N.vals <- c(50000, 100000, 200000, 300000, 400000)
pop <- simulate.population(M = 10000, h2 = 0.5, n.target)
var.vs.N <- sweep.var(values = N.vals, function(N) prs.posterior.var(pop, N, K)$pgs.var)
plot.N <- plot.violin(var.vs.N, "N", add.sep) #, expected = function(N) expected.post.var(M = 10000, h2 = 0.5, N = N, K = K))
print(plot.N)

ggsave(here("figures", "var-vs-n.png"), plot = plot.N, width = 8, height = 8, dpi = 300)


# Effect of M on prs var: independent population for each M
set.seed(111)
#M.vals <- c(100, 1000, 10000, 100000)
M.vals <- c(1000, 10000, 50000, 100000)
var.vs.M <- sweep.var(M.vals, function(M) prs.posterior.var(simulate.population(M, h2 = 0.5, n.target), 100000, K)$pgs.var)
#plot.M <- plot.violin(var.vs.M, "M", add.sep, expected = function(M) expected.post.var(M = M, h2 = 0.5, N = 100000, K = K), log.y = TRUE)
plot.M <- plot.violin(var.vs.M, "M", add.sep) #, expected = function(M) expected.post.var(M = M, h2 = 0.5, N = 100000, K = K))
print(plot.M)

ggsave(here("figures", "var-vs-n.png"), plot = plot.M, width = 8, height = 8, dpi = 300)


# Effect of h2: independent population for each h2
set.seed(111)
h2.vals <- c(0.1, 0.2, 0.4, 0.6, 0.8)
#h2.vals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
var.vs.h2 <- sweep.var(h2.vals, function(h2) prs.posterior.var(simulate.population(10000, h2, n.target), 100000, K)$pgs.var)
plot.h2 <- plot.violin(var.vs.h2, expression(h^2)) #, expected = function(h2) expected.post.var(M = 10000, h2 = h2, N = 100000, K = K))
print(plot.h2)

ggsave(here("figures", "var-vs-h2.png"), plot = plot.h2, width = 8, height = 8, dpi = 300)
