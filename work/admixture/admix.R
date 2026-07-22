library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Sample data structure (assuming you have data like this)
# Your data should have columns: Sample, Group, K, and ancestry components (e.g., Pop1, Pop2, Pop3, Pop4, Pop5)

# Function to prepare admixture data with consistent ordering
prepare_admixture_data <- function(data_k4, data_k5) {
  # Step 1: Determine ordering based on K=4 most dominant component
  order_df <- data_k5 %>%
    # Find dominant component for each sample in K=4
    pivot_longer(
      cols = starts_with("K"),
      names_to = "Component",
      values_to = "Proportion"
    ) %>%
    group_by(samples) %>%
    slice_max(Proportion, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Create ordering within each group
    arrange(previous_molecular_id, Component, desc(Proportion)) %>%
    mutate(Order = row_number()) %>%
    select(samples, Order)

  # Step 2: Apply this ordering to both datasets
  data_k4_ordered <- data_k4 %>%
    left_join(order_df, by = "samples") %>%
    arrange(Order) %>%
    mutate(samples = factor(samples, levels = unique(samples)))

  data_k5_ordered <- data_k5 %>%
    left_join(order_df, by = "samples") %>%
    arrange(Order) %>%
    mutate(samples = factor(samples, levels = levels(data_k4_ordered$samples)))

  return(list(k4 = data_k4_ordered, k5 = data_k5_ordered))
}

# Function to reshape data for plotting
reshape_for_plotting <- function(data, k_value) {
  data %>%
    pivot_longer(
      cols = starts_with("K"),
      names_to = "Ancestry",
      values_to = "Proportion"
    ) %>%
    mutate(K = as.factor(k_value)) |>
    filter(Ancestry != "key") |>
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
}

# Example data preparation (replace with your actual data)
# Assuming your data has columns: Sample, Group, Pop1, Pop2, Pop3, Pop4 (and Pop5 for K=5)
k4_props <- read_table(
  "K4/LD_pruned.MN10.MX100.MS80.MQ20.admixture_input.4.19.Q",
  col_names = c("K1", "K2", "K3", "K4")
)
k5_props <- read_table(
  "K5/LD_pruned.MN10.MX100.MS80.MQ20.admixture_input.5.16.Q",
  col_names = c("K1", "K2", "K3", "K4", "K5")
)
meta <- read_csv("../../data/samples.csv")
samples <- read_lines("../samples_selected.txt")

data_k4 <- k4_props |>
  cbind(samples) |>
  left_join(meta, by = c("samples" = "new_id")) |>
  distinct(samples, .keep_all = TRUE)

data_k5 <- k5_props |>
  cbind(samples) |>
  left_join(meta, by = c("samples" = "new_id")) |>
  distinct(samples, .keep_all = TRUE)

# Prepare ordered data
ordered_data <- prepare_admixture_data(data_k4, data_k5)

# Reshape for plotting
plot_data_k4 <- reshape_for_plotting(ordered_data$k4, 4)
plot_data_k5 <- reshape_for_plotting(ordered_data$k5, 5)

plot_data_k4 <- plot_data_k4 |>
  mutate(Ancestry = factor(Ancestry, c("K1", "K2", "K4", "K3")))

# Alternative: If you want separate plots with consistent ordering
p4 <- ggplot(plot_data_k4, aes(x = samples, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 1, show.legend = FALSE) +
  facet_wrap(
    ~previous_id,
    scales = "free_x",
    space = "free_x",
    strip.position = "bottom"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "K = 4") +
  scale_fill_manual(
    values = c(
      "#e41a1c",
      "#4daf4a",
      "#984ea3",
      "#ff7f00"
    )
  ) +
  scale_y_continuous(expand = c(0, 0))

p5 <- ggplot(plot_data_k5, aes(x = samples, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 1, show.legend = FALSE) +
  facet_wrap(
    ~previous_id,
    scales = "free_x",
    space = "free_x",
    strip.position = "bottom"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = "K = 5") +
  scale_fill_manual(
    values = c(
      "#e41a1c",
      "#377eb8",
      "#984ea3",
      "#4daf4a",
      "#ff7f00"
    )
  ) +
  scale_y_continuous(expand = c(0, 0))

# Display plots
library(patchwork)
composite_plot <- p4 / p5
ggsave("Fig3_admixture.png", composite_plot, width = 9, height = 6)
