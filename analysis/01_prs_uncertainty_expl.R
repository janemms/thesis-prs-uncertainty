library(here)
source(here("R", "thesis-functions.R"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ----------------- scatter plot of PRS vs PRS SD and --------------------------
# ---------- PRS densities of affected and unaffected individuals --------------

set.seed(111)

pop <- simulate.population(M = 10000, h2 = 0.5, n.target = 1000)
case.status <- simulate.case.status(pop, h2 = 0.5, K = 0.2)
prs <- prs.posterior.var(pop, N = 100000, K = 0.2)

plot.df <- data.frame(prs = prs$pgs, prs.var = prs$pgs.var, case.status = factor(case.status, labels = c("Unaffected", "Affected")))

# compute metrics
var.summary <- c(
  mean = mean(plot.df$prs.var, na.rm = TRUE),
  sd = sd(plot.df$prs.var, na.rm = TRUE),
  CV = sd(plot.df$prs.var, na.rm = TRUE)/mean(plot.df$prs.var, na.rm = TRUE),
  min = min(plot.df$prs.var, na.rm = TRUE),
  max = max(plot.df$prs.var, na.rm = TRUE),
  median = median(plot.df$prs.var, na.rm = TRUE)
)
print(var.summary)

cor.prs.var <- cor.test(plot.df$prs, plot.df$prs.var, method = "pearson")
print(cor.prs.var)

cor.disease.var <- cor.test(as.integer(plot.df$case.status == "Affected"), plot.df$prs.var, method = "pearson")
print(cor.disease.var)

# Summary statistics by disease status
aggregate(
  prs ~ case.status,
  data = plot.df,
  FUN = function(x) c(
    mean = mean(x),
    sd = sd(x),
    median = median(x)
  )
)

# Difference in mean PRS and 95% confidence interval
prs.mean.test <- t.test(prs ~ case.status, data = plot.df)
print(prs.mean.test)

# Standardized mean difference (Cohen's d)
prs.groups <- split(plot.df$prs, plot.df$case.status)

mean.diff <- mean(prs.groups$Affected) -
  mean(prs.groups$Unaffected)

pooled.sd <- sqrt(
  ((length(prs.groups$Affected) - 1) * var(prs.groups$Affected) +
     (length(prs.groups$Unaffected) - 1) * var(prs.groups$Unaffected)) /
    (length(prs.groups$Affected) + length(prs.groups$Unaffected) - 2)
)

cohens.d <- mean.diff / pooled.sd

cat("Mean difference (affected - unaffected):", mean.diff, "\n")
cat("Cohen's d:", cohens.d, "\n")

# plot

disease.colors <- c("Unaffected" = "#9ECAE1", "Affected" = "#003999")

prs.vs.var <- ggplot(plot.df, aes(x=prs, y=prs.var, color=case.status)) +
  geom_point(alpha = 0.6, size=1.2) +
  theme_bw() +
  scale_color_manual(values = disease.colors) +
  guides(color = "none") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.border = element_blank()
  ) +
  labs(x = "PRS", y = "Posterior variance", color = "Disease status")
print(prs.vs.var)
ggsave(filename = here("figures", "prs-vs-var.png"), plot = prs.vs.var, width = 8, height = 8, dpi = 300)

aff.vs.unaff <- ggplot(plot.df, aes(x=prs, fill= case.status)) +
  geom_density(alpha = 0.6, size = 0.5) +
  scale_fill_manual(values = disease.colors) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.border = element_blank()
  ) +
  labs(x = "PRS ", y= "Density", fill = "Disease status") 
print(aff.vs.unaff)

ggsave(here("figures", "affected-vs-unaffected.png"), plot = aff.vs.unaff, width = 8, height = 8, dpi = 300)

expl.combined <- prs.vs.var + aff.vs.unaff + plot_layout(guides = "collect") & theme(
  legend.position = "bottom",
  legend.direction = "horizontal", # Optional: lays the legend items out horizontally
  legend.title = element_text(face = "bold", size = 11),
  legend.text = element_text(size = 10)
)
print(expl.combined)
ggsave(here("figures", "explorative-combined.png"), plot = expl.combined, width = 8, height = 5, dpi = 300)

