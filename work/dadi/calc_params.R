library(tidyverse)

data <- read_csv("best_model_comparison.csv")

mutation <- 7.9e-9

L <- (234.9 * 1107077) * (1065048 / 206788238)

data <- data |>
  mutate(
    nref = theta / (4 * mutation * L),
    real_nu1 = nref * nu1,
    real_nu2 = nref * nu2,
    real_m = m / (2 * nref),
    real_m12 = m12 / (2 * nref),
    real_m21 = m21 / (2 * nref),
    real_t1 = T1 * 2 * nref,
    real_t2 = T2 * 2 * nref,
    tsplit = real_t1 + real_t2
  )
write_csv(data, "real_params.csv")

data |>
  filter(model == "asym_mig") |>
  select(comparison, real_t1, real_m12, real_m21) |>
  pivot_longer(real_m12:real_m21) |>
  group_by(comparison) |>
  summarise(real_t1 = first(real_t1), value = max(value)) |>
  ggplot(aes(x = real_t1, y = value, colour = comparison)) +
  geom_point(size = 3) +
  scale_colour_brewer(palette = "Paired", name = "Comparison") +
  theme_minimal()

asym_mig <- data |>
  filter(model == "asym_mig") |>
  ggplot(
    aes(x = real_t1, y = real_m12 + real_m21, colour = comparison)
  ) +
  geom_point(size = 3) +
  labs(
    x = "Time of split (generations)",
    y = "Total migration rate (m12 + m21)"
  ) +
  scale_colour_brewer(palette = "Paired", name = "Comparison") +
  theme_minimal()
ggsave("asym_mig.divergence_vs_migration.png", asym_mig)

sym_mig <- data |>
  filter(model == "split_mig") |>
  ggplot(
    aes(x = real_t1, y = real_m, colour = comparison)
  ) +
  geom_point(size = 3) +
  labs(
    x = "Time of split (generations)",
    y = "Migration rate"
  ) +
  scale_colour_brewer(palette = "Paired", name = "Comparison") +
  theme_minimal()
ggsave("split_mig.divergence_vs_migration.png", sym_mig)
