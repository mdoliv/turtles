library(tidyverse)

refs <- c(
  "rCarCar",
  "rCheMyd",
  "rDerCor",
  "rEreImb",
  "rLepKem",
  "rLepOli",
  "rNatDep"
)

full_data <- data.frame()

for (ref in refs) {
  ref_data <- read_tsv(
    paste0(ref, ".filtered.depths.tsv"),
    col_names = c(
      "sample",
      "n_loci",
      "n_used_fw_reads",
      "mean_cov",
      "mean_cov_ns"
    )
  )
  ref_data$ref <- ref
  full_data <- rbind(full_data, ref_data)
}

#meta <- read_tsv("sample_pops.tsv", col_names = c("sample", "pop"))
#
#dcor_data <- full_data |>
#  filter(ref == "rDerCor") |>
#  separate_wider_delim(
#    sample,
#    delim = ".",
#    names = c("sample", "method", "lib", "mod")
#  )
#
#meta |>
#  left_join(dcor_data, by = "sample") |>
#  group_by(pop) |>
#  slice_max(mean_cov_ns, n = 2)

#full_data |>
#  summary()

meta <- read_csv("../../data/samples.csv") |>
  distinct(new_id, .keep_all = TRUE)

full_data <- full_data |>
  mutate(sample = str_remove(sample, "\\..*")) |>
  left_join(meta, by = c("sample" = "new_id"))

coverage_plot <- full_data |>
  ggplot(aes(
    x = ref,
    y = mean_cov_ns,
    fill = previous_molecular_id,
    colour = previous_molecular_id
  )) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Reference", y = "Mean weighted depth") +
  scale_fill_manual(
    name = "Previous molecular\nidentification",
    labels = c("Hybrid", "Cc", "Cm", "Ei", "Lo", "Non identified"),
    values = c("darkgrey", "#F28E2B", "#59A14F", "#4E79A7", "#B07AA1"),
    na.value = "black"
  ) +
  scale_colour_manual(
    name = "Previous molecular\nidentification",
    labels = c("Hybrid", "Cc", "Cm", "Ei", "Lo", "Non identified"),
    values = c("darkgrey", "#F28E2B", "#59A14F", "#4E79A7", "#B07AA1"),
    na.value = "black"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())
#coverage_plot <- full_data |>
#  ggplot(aes(
#    x = fct_reorder(sample, mean_cov_ns),
#    y = ref,
#    fill = mean_cov_ns
#  )) +
#  geom_tile() +
#  scale_fill_viridis_c(
#    name = "Cobertura média\nponderada",
#    option = "magma"
#  ) +
#  labs(x = "Amostra", y = NULL) +
#  theme_minimal() +
#  theme(panel.grid = element_blank(), axis.text.x = element_blank())

ggsave("coverage_plot.png", coverage_plot, width = 8, height = 6)
