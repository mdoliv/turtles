# =============================================================================
# ADMIXTURE Results Plotter
# Plots K=4 and K=5 results with samples ordered within each species by their
# proportion of the dominant K=5 component for that species.
# =============================================================================

library(tidyverse)
library(patchwork)

# =============================================================================
# 1. LOAD YOUR DATA — edit these paths and column names
# =============================================================================

# --- Q matrices (one row per sample) ---
# Replace with your actual file paths. Common formats:
#   read.table("file.5.Q")           # plain ADMIXTURE output
#   read_csv("file_k5.csv")          # CSV with or without header

q4_raw <- read.table("K4/LD_pruned.MN10.MX100.MS80.MQ20.admixture_input.4.19.Q") # K=4 Q matrix
q5_raw <- read.table("K5/LD_pruned.MN10.MX100.MS80.MQ20.admixture_input.5.16.Q") # K=5 Q matrix

# --- Sample metadata ---
# Must contain at least: sample ID and species columns.
# Rows must be in the same order as the Q matrices.
meta <- read_csv("../../data/samples.csv")
kept_samples <- read_lines("../samples_selected.txt")
meta <- meta |>
  filter(new_id %in% kept_samples) |>
  distinct(new_id, .keep_all = TRUE)
# Expected columns (rename as needed below):
#   sample_id : unique identifier per sample
#   species   : grouping variable for ordering and faceting

# Rename to standard names if your columns differ:
meta <- meta %>%
  rename(
    sample_id = new_id, # e.g. "Sample", "ID", "ind"
    species = previous_molecular_id # e.g. "Species", "Population", "pop"
  ) |>
  mutate(
    species = case_when(
      species == "admixed" ~ "Híbrido",
      species == "cc" ~ "Cc",
      species == "cm" ~ "Cm",
      species == "ei" ~ "Ei",
      species == "lo" ~ "Lo",
      is.na(species) ~ "Não identificado",
      .default = species
    ),
    species = factor(
      species,
      levels = c("Cc", "Ei", "Cm", "Lo", "Híbrido", "Não identificado")
    )
  )

# =============================================================================
# 2. PREPARE LONG-FORMAT DATA FRAMES
# =============================================================================

prepare_q <- function(q_raw, meta, k_label) {
  colnames(q_raw) <- paste0("C", seq_len(ncol(q_raw)))
  q_raw %>%
    bind_cols(meta %>% select(sample_id, species)) %>%
    pivot_longer(
      cols = starts_with("C"),
      names_to = "component",
      values_to = "proportion"
    ) %>%
    mutate(k = k_label)
}

q5_long <- prepare_q(q5_raw, meta, "K = 5")
q4_long <- prepare_q(q4_raw, meta, "K = 4")

# =============================================================================
# 3. DETERMINE SAMPLE ORDER FROM K=5
# =============================================================================

# Step 1 – For each species, find which K=5 component is most dominant
#          (highest mean proportion across all samples in that species).
dominant_component_per_species <- q5_long %>%
  group_by(species, component) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop") %>%
  slice_max(mean_prop, by = species, n = 1) %>%
  select(species, dominant_component = component)

# Step 2 – Within each species, sort samples by their proportion of that
#          dominant component (descending).
sample_order <- q5_long %>%
  left_join(dominant_component_per_species, by = "species") %>%
  filter(component == dominant_component) %>%
  arrange(species, desc(proportion)) %>%
  pull(sample_id)

# Apply the order as a factor level — both K panels will respect this.
q5_long <- q5_long %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order),
    component = factor(component, levels = c("C1", "C2", "C3", "C5", "C4"))
  )
q4_long <- q4_long %>%
  mutate(
    sample_id = factor(sample_id, levels = sample_order),
    component = factor(component, levels = c("C1", "C2", "C4", "C3"))
  )

# =============================================================================
# 4. COLOUR PALETTE
# =============================================================================

# Earthy, distinct palette that reads well in print and on screen.
# Edit freely — must have at least max(K4, K5) colours.
component_colours <- c(
  C1 = "#4E79A7", # steel blue
  C2 = "#F28E2B", # amber
  C3 = "#59A14F", # moss green
  C4 = "#E15759", # coral red
  C5 = "#B07AA1" # muted violet
)

# =============================================================================
# 5. PLOT FUNCTION
# =============================================================================

p4 <- ggplot(q4_long, aes(x = sample_id, y = proportion, fill = component)) +
  geom_col(width = 1, position = "stack", linewidth = 0, show.legend = FALSE) +
  facet_grid(
    cols = vars(species),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "#4E79A7",
      "#59A14F",
      "#B07AA1",
      "#F28E2B"
    )
  ) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format()) +
  labs(title = "K = 4", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    # --- axes ---
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "grey30"),
    axis.title.y = element_text(
      size = 9,
      color = "grey30",
      margin = margin(r = 6)
    ),
    # --- facet strips ---
    strip.text = element_text(size = 9, color = "grey20"),
    strip.background = element_blank(),
    # --- panel ---
    panel.grid = element_blank(),
    panel.spacing.x = unit(2, "pt"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    # --- legend ---
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(size = 8),
    # --- titles ---
    plot.title = element_text(
      size = 11,
      color = "grey15",
      margin = margin(b = 4)
    ),
    plot.margin = margin(4, 6, 4, 6)
  )

p5 <- ggplot(q5_long, aes(x = sample_id, y = proportion, fill = component)) +
  geom_col(width = 1, position = "stack", linewidth = 0, show.legend = FALSE) +
  facet_grid(
    cols = vars(species),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      "#4E79A7",
      "#E15759",
      "#B07AA1",
      "#F28E2B",
      "#59A14F"
    )
  ) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format()) +
  labs(title = "K = 5", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    # --- axes ---
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "grey30"),
    axis.title.y = element_text(
      size = 9,
      color = "grey30",
      margin = margin(r = 6)
    ),
    # --- facet strips ---
    strip.text = element_text(size = 9, color = "grey20"),
    strip.background = element_blank(),
    # --- panel ---
    panel.grid = element_blank(),
    panel.spacing.x = unit(2, "pt"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    # --- legend ---
    legend.position = "right",
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(size = 8),
    # --- titles ---
    plot.title = element_text(
      size = 11,
      color = "grey15",
      margin = margin(b = 4)
    ),
    plot.margin = margin(4, 6, 4, 6)
  )

# =============================================================================
# 6. COMBINE AND SAVE
# =============================================================================

combined <- p4 /
  p5 +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = "admixture_k4_k5.pdf",
  plot = combined,
  width = 9,
  height = 6, # adjust to taste
  device = cairo_pdf
)

# Also save as PNG for quick preview
ggsave(
  filename = "admixture_k4_k5.png",
  plot = combined,
  width = 9,
  height = 6,
  dpi = 300
)

message("Plots saved: admixture_k4_k5.pdf / .png")
