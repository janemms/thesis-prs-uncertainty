library(here)
source(here("R", "thesis-functions.R"))

library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(111)

M <- 10000
n.target <- 10000
nreps <- 10
K = 0.1
t <- qnorm(1-K) # threshold for disease

# parameter combinations
test.grid <- expand.grid(N = c(50000, 200000), h2 = c(0.2, 0.5, 0.8))
#summary.list <- vector("list", nrow(test.grid))

# one rep: single population and single gwas draw
run.rep <- function(row, rep) {
  N <- test.grid$N[row]
  h2 <- test.grid$h2[row]
  set.seed(get.seed(row, rep))  # 
  pop <- simulate.population(M, h2, n.target)
  prs <- prs.posterior.var(pop, N, K)
  
  data.frame(N, h2, btw = var(prs$pgs), within = mean(prs$pgs.var))
}

all.runs <- expand.grid(row = 1:nrow(test.grid), rep = 1:nreps) # table with all combinations of row and rep
reps <- do.call(rbind, Map(run.rep, all.runs$row, all.runs$rep)) # call run.rep once for each element in all.runs

# avg over reps
summary.df <- reps %>%
  group_by(N, h2) %>%
  summarise(btw = mean(btw), within = mean(within), .groups = "drop") # remove all grouping

print(mutate(summary.df, total = btw + within))

# plot
decomp.df <- summary.df %>%
  pivot_longer(cols = c(btw, within),
               names_to = "component", values_to = "value") %>%
  group_by(N, h2) %>%
  mutate(prop = value / sum(value)) %>%
  ungroup() %>%
  mutate(component = factor(component, levels = c("btw", "within"),
                            labels = c("Between-individual",
                                       "Within-individual")),
                            N = factor(N))

plot.var.decomp <- ggplot(decomp.df, aes(x = N, y = value, fill = component)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.8, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * prop)), position = position_stack(vjust = 0.5, reverse = TRUE), color = "white", size = 3.5, fontface = "bold") +
  facet_wrap(~h2, labeller = as_labeller(function(x) paste0("h\u00B2 = ", x))) +
  #scale_color_manual(values = comp.cols) +
  scale_fill_manual(values = c("Between-individual" = "#56B4E9", "Within-individual" = "#0072B2")) +
  scale_x_discrete(labels = function(x) formatC(as.numeric(x), format = "d", big.mark = ","))+
  labs(x = "GWAS sample size", y = "Variance", fill = "Component") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5), 
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        strip.text         = element_text(face = "bold", size = 13),
        legend.title       = element_text(face = "bold", size = 13),
        legend.text         = element_text(size = 13)
  ) 

print(plot.var.decomp)
ggsave(here("figures", "var-decomposition.png"), plot = plot.var.decomp, width = 8, height = 5, dpi = 300)
