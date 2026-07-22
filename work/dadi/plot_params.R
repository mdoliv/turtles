library(tidyverse)

params_ci <- read_csv("params_with_ci.csv")

ci_plot <- params_ci |>
  ggplot(aes(x = model, y = value, colour = comparison)) +
  geom_pointrange(
    aes(ymax = upper, ymin = lower),
    position = position_dodge(width = 1),
    size = 0.8
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  labs(x = NULL, y = NULL) +
  facet_wrap(~param, scales = "free_y") +
  scale_color_brewer(palette = "Paired", name = NULL) +
  scale_x_discrete(
    labels = c("Asymmetric\nmigration", "Symmetric\nmigration")
  ) +
  scale_y_continuous(labels = scales::label_number()) +
  theme_minimal() +
  theme(
    strip.background = element_rect(colour = "black"),
    legend.position = c(0.85, 0.3),
    legend.background = element_rect(colour = "black")
  )
ggsave("params_ci.png", ci_plot, width = 12, height = 9)
