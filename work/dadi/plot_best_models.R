library(tidyverse)

best_models <- read_csv("best_model_comparison.csv")
best_models <- best_models |>
  mutate(
    model = factor(
      model,
      c(
        "no_mig",
        "split_mig",
        "sec_contact_sym_mig",
        "asym_mig",
        "sec_contact_asym_mig"
      )
    ),
    significant = GIM_LRT_p < 0.05
  )

model_plot <- best_models |>
  ggplot(aes(
    x = model,
    y = likelihood,
    colour = comparison,
    group = comparison
  )) +
  geom_line() +
  geom_point(size = 2, aes(shape = significant)) +
  guides(shape = "none") +
  scale_x_discrete(
    labels = c(
      "No migration",
      "Symmetric\nmigration",
      "Secondary\ncontact",
      "Asymmetric\nmigration",
      "Sec. contact\nasym. migration"
    )
  ) +
  labs(x = NULL, y = "Log likelihood") +
  scale_colour_brewer(palette = "Paired", name = NULL) +
  scale_shape_manual(values = c(1, 16), na.value = 1) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "top")

ggsave("best_models.png", model_plot)
