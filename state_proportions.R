library(tidyverse)

sample_meta <- read_csv("data/samples.csv")
kept_samples <- read_lines("data/samples_depth5.txt")
sample_meta <- sample_meta |>
  filter(new_id %in% kept_samples) |>
  distinct(new_id, .keep_all = TRUE) |>
  mutate(
    type = case_when(
      previous_molecular_id %in% c("cc", "ei", "cm", "lo") ~ "Pure",
      previous_molecular_id == "admixed" ~ "Hybrid",
      .default = "Not identified"
    ),
    state = case_when(
      state == "Southeast Atlantic" ~ "Elevação do Rio Grande",
      .default = state
    )
  )

proportion_plot <- sample_meta |>
  ggplot(aes(y = state, fill = type)) +
  geom_bar(position = "fill", colour = "black") +
  scale_fill_manual(
    values = c("darkgrey", "black", "white"),
    name = "Previous molecular\nidentification"
  ) +
  scale_x_continuous(expand = FALSE, labels = scales::label_percent()) +
  labs(x = "Proportion", y = "State") +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggsave("sample_proportions.png", proportion_plot)
