library(tidyverse)

runs_summary <- read_csv("summary.csv")
runs_summary <- runs_summary |>
  group_by(k) |>
  slice_min(cv) |>
  slice_max(ll)

best_k_plot <- runs_summary |>
  ggplot(aes(x = as.factor(k), y = cv)) +
  geom_vline(
    xintercept = runs_summary[which.min(runs_summary$cv), ]$k,
    linetype = "dashed"
  ) +
  geom_line(group = 1) +
  geom_point(size = 3, colour = "red") +
  labs(y = "CV error", x = "K (best run)") +
  theme_minimal()
ggsave("best_k.png", best_k_plot)

best_run <- read_table(
  "K5/LD_pruned.MN10.MX100.MS80.MQ20.admixture_input.5.16.Q",
  col_names = FALSE
)
samples_list <- read_table("../pixy/popmap.tsv", col_names = FALSE)
meta <- read_csv("../../data/samples.csv")

run_with_samples <- cbind(best_run, samples_list)
names(run_with_samples) <- c("K1", "K2", "K3", "K4", "K5", "sample", "pop")

run_with_samples <- run_with_samples |>
  left_join(meta, by = c("sample" = "new_id")) |>
  distinct(sample, .keep_all = TRUE)

to_plot <- run_with_samples |>
  pivot_longer(
    K1:K5,
    names_to = "comp",
    values_to = "prop"
  ) |>
  mutate(
    previous_id = case_when(
      is.na(previous_molecular_id) ~ "Not identified",
      previous_molecular_id == "admixed" ~ "Hybrid/admixed",
      previous_molecular_id == "cc" ~ "Cc",
      previous_molecular_id == "cm" ~ "Cm",
      previous_molecular_id == "ei" ~ "Ei",
      previous_molecular_id == "lo" ~ "Lo"
    ),
    previous_id = factor(
      previous_id,
      c("Cc", "Ei", "Cm", "Lo", "Hybrid/admixed", "Not identified")
    )
  )

dominant_component <- to_plot |>
  group_by(previous_id, comp) |>
  summarise(mean_prop = mean(prop), .groupds = "drop") |>
  group_by(previous_id) |>
  slice_max(mean_prop, n = 1) |>
  select(previous_id, dominant_component = comp)

to_plot <- to_plot |>
  left_join(dominant_component, by = c("previous_id"))

dominant_props <- to_plot |>
  filter(comp == dominant_component) |>
  select(previous_id, sample, dominant_prop = prop)

to_plot <- to_plot |>
  left_join(dominant_props, by = c("previous_id", "sample"))

to_plot <- to_plot |>
  group_by(previous_id) |>
  mutate(sample_ordered = reorder(sample, -dominant_prop)) |>
  ungroup()

admixture_plot <- to_plot |>
  ggplot(
    aes(x = sample_ordered, y = prop, fill = comp)
  ) +
  geom_col(show.legend = FALSE, width = 1) +
  facet_wrap(
    ~previous_id,
    scales = "free_x",
    space = "free_x",
    strip.position = "bottom"
  ) +
  theme_minimal() +
  labs(y = NULL) +
  scale_y_continuous(expand = NA) +
  scale_fill_manual(
    values = c(
      "#e41a1c",
      "#377eb8",
      "#984ea3",
      "#4daf4a",
      "#ff7f00"
    )
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.05, "cm"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.title.x = element_blank()
  )
ggsave("Fig3_admixture.png", admixture_plot, width = 8, height = 4)

###

best_run_k4 <- read_table(
  "K4/LD_pruned.MN10.MX100.MS80.MQ20.admixture_input.4.19.Q",
  col_names = FALSE
)
k4_with_samples <- cbind(best_run_k4, samples_list)
names(k4_with_samples) <- c("K1", "K2", "K3", "K4", "sample", "pop")

k4_with_samples <- to_plot |>
  left_join(k4_with_samples, by = c("sample_ordered" = "sample")) |>
  select(-c(prop, comp))

to_plot_k4 <- k4_with_samples |>
  pivot_longer(
    K1:K4,
    names_to = "comp",
    values_to = "prop"
  ) |>
  mutate(
    comp = factor(
      comp,
      c("K1", "K2", "K4", "K3")
    ),
    sample = factor(
      sample,
      levels = levels(to_plot$sample_ordered)
    )
  )

k4_plot <- to_plot_k4 |>
  ggplot(
    aes(x = sample_ordered, y = prop, fill = comp)
  ) +
  geom_col(show.legend = FALSE, width = 1) +
  facet_wrap(
    ~previous_id,
    scales = "free_x",
    space = "free_x",
    strip.position = "bottom"
  ) +
  theme_minimal() +
  labs(y = NULL) +
  scale_y_continuous(expand = NA) +
  scale_fill_manual(
    values = c(
      "#e41a1c",
      "#4daf4a",
      "#984ea3",
      "#ff7f00"
    )
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.05, "cm"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.title.x = element_blank()
  )
ggsave("Fig3_admixture.png", admixture_plot, width = 8, height = 4)

library(patchwork)
k4_plot / admixture_plot


###
species_delim <- run_with_samples |>
  mutate(
    spp_assignment = case_when(
      K1 >= 0.95 | K4 >= 0.95 ~ "cc",
      K2 >= 0.95 ~ "ei",
      K3 >= 0.95 ~ "lo",
      K5 >= 0.95 ~ "cm",
      .default = "admixed"
    )
  ) |>
  select(sample, spp_assignment)

write_csv(species_delim, "species_delim.csv")

run_with_samples |>
  select(K1:sample) |>
  write_csv("samples_with_k.csv")
