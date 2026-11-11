
library(ggplot2)

# Define distribution
mu <- 0
sigma <- 1
t <- 0.5

# Generate data
df <- data.frame(x = seq(-4, 4, length.out = 1000))
df$y <- dnorm(df$x, mean = mu, sd = sigma)
df$shade <- ifelse(df$x > t, "Above threshold", "Below threshold")

# Plot
ggplot(df, aes(x, y)) +
  geom_area(aes(fill = shade), alpha = 0.4) +
  geom_line(linewidth = 1.1, color = "black") +
  geom_vline(xintercept = t, linetype = "dashed",
             linewidth = 0.8) +
  annotate("text", x = t + 1, y = 0.4,
           label = "Threshold t") +
  scale_fill_manual(values = c(
    "Above threshold" = "#4C9F70",
    "Below threshold" = "white"
  )) +
  labs(
    x = "Genetic value (GV)",
    y = "Density",
    fill = ""
  ) +
  ylim(0, NA) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

