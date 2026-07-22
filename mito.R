library(tidyverse)

sample_meta <- read_csv("data/samples.csv")
kept_samples <- read_lines("data/samples_depth5.txt")
sample_meta <- sample_meta |>
  filter(new_id %in% kept_samples) |>
  distinct(new_id, .keep_all = TRUE)

mito <- sample_meta |>
  filter_out(is.na(mtDNA)) |>
  filter(previous_molecular_id == "admixed") |>
  select(new_id, crossing, mtDNA)

print(mito)
