library(ggplot2)
library(dplyr)


# run one simulation per parameter value
sweep.var <- function(values, one.run) {
  mat <- sapply(values, one.run)
  colnames(mat) <- values
  as.data.frame(mat) %>%
    pivot_longer(everything(), names_to = "value", values_to = "post.var") %>%
    mutate(value = factor(value, levels = as.character(values)))
}


# create violin plots
plot.violin <- function(df, x.label, x.labels = waiver()){ # waiver: works for h2 too
  p <- ggplot(df, aes(x=value, y=post.var)) +
    geom_violin(alpha = 0.35, trim = FALSE, fill = "#0072B2", colour = "#0072B2") +
    geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.6, fill = "#0072B2", colour = "#0072B2") +
    scale_x_discrete(labels = x.labels) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_bw() + 
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.8),
          panel.border = element_blank()
    ) +
    labs(x = x.label, y = "Posterior variance")
  p
}

# create scatter plots
plot.risk.scatter <- function(df) {
  ggplot(df, aes(x = p.m1, y = p.m2)) +
    geom_point(alpha = 0.35, size = 0.7, color = "#0072B2") +
    geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
    facet_grid(cols = vars(panel)) +
    theme_classic() +
    theme(
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(x = "Disease probability - Method P", y = "Disease probability - Method U")
}


# plot method comparisons as boxplots
pastel.colors <- c(
  "Method P" = "#9ECAE1",
  "Method U" = "#0072B2"
)

scenarios.keep <- data.frame(
  cost.FN = c(10, 3, 1, 1, 1),
  cost.FP = c(1, 1, 1, 3, 10)
)

prepare.results <- function(df, M.val, K.val, N.val, h2.val) {
  df %>% 
    filter(M %in% M.val, N %in% N.val, K %in% K.val, h2 %in% h2.val, 
           method %in% c("point.estimate", "uncertainty.aware")) %>%
    semi_join(scenarios.keep, by = c("cost.FN", "cost.FP")) %>%
    mutate(method = recode(method,
                           "point.estimate" = "Method P",
                           "uncertainty.aware" = "Method U"),
           prop.selected = n.selected / n, 
           cost.per.ind = cost / n,
           cost.ratio = cost.FN / cost.FP,
           scenario = factor(paste0(cost.FN, ":", cost.FP),
                             levels = unique(paste0(cost.FN, ":", cost.FP))[order(unique(cost.ratio))]))
}

plot.comparison.boxplot <- function(df, facet.var = c("N", "h2"), y.var, y.label){
  facet.var <- match.arg(facet.var)
  facet.label <- switch(facet.var, 
                        N = function(x) paste0("N = ", format(as.numeric(x), big.mark = ",", trim = TRUE, scientific = FALSE)),
                        h2 = function(x) paste0("h\u00B2 = ", x))
  
  ggplot(df, aes(x = scenario, y = .data[[y.var]], fill = method, color = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_grid(cols = vars(.data[[facet.var]]), labeller = as_labeller(facet.label)) +
    scale_color_manual(values = pastel.colors, name = "Method") +
    scale_fill_manual(values = pastel.colors, name = "Method") +
    theme_classic() +
    theme(panel.border = element_blank(),        
          panel.grid.major.x = element_blank(),  
          panel.grid.minor.y = element_line(color = "grey85", linewidth = 0.5),
          panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5), 
          axis.line = element_line(color = "black", linewidth = 0.5), 
          strip.text         = element_text(face = "bold", size = 12),
          legend.title       = element_blank(),#element_text(face = "bold", size = 11),
          legend.text         = element_text(size = 10),
          legend.position = "bottom") +
    expand_limits(y = 0) +
    labs(x = expression("Cost scenario (" * C[FN]:C[FP]*")"), y = y.label)
}