library(rnaturalearth)
library(geobr)
library(tidyverse)
library(ggtext)
library(terra)
library(scatterpie)

sample_meta <- read_csv("data/samples.csv")
kept_samples <- read_lines("data/samples_depth5.txt")
sample_meta <- sample_meta |>
  filter(new_id %in% kept_samples) |>
  distinct(new_id, .keep_all = TRUE) |>
  mutate(
    lat = case_when(
      is.na(capture_latitude) ~ imputed_latitude,
      .default = capture_latitude
    ),
    lon = case_when(
      is.na(capture_longitude) ~ imputed_longitude,
      .default = capture_longitude
    )
  )

scatterpie_data <- sample_meta |>
  mutate(
    type = case_when(
      previous_molecular_id %in% c("cc", "ei", "cm", "lo") ~ "Pure",
      previous_molecular_id == "admixed" ~ "Hybrid",
      .default = "Not identified"
    )
  ) |>
  group_by(type, lat, lon) |>
  summarise(count = n()) |>
  pivot_wider(names_from = type, values_from = count, values_fill = 0)

ne1 <- terra::rast("NE1_50M_SR_W/NE1_50M_SR_W.tif")
bbox <- terra::ext(-80, -25, -40, 25)
ne1_crop <- terra::crop(ne1, bbox)

countries <- ne_countries(type = "countries", scale = "medium")
states <- ne_states(country = "Brazil")
countries_sa <- sf::st_crop(
  countries,
  xmin = -80,
  xmax = -25,
  ymin = -40,
  ymax = 25
)

countries_sa |>
  ggplot() +
  tidyterra::geom_spatraster_rgb(data = ne1_crop) +
  geom_sf(fill = NA) +
  geom_sf(data = states, fill = NA) +
  geom_scatterpie(
    data = scatterpie_data,
    mapping = aes(x = lon, y = lat),
    cols = c("Hybrid", "Not identified", "Pure")
  ) +
  coord_sf(xlim = c(-80, -25), ylim = c(-40, 25), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    legend.text = element_markdown(size = 12),
    legend.title = element_text(size = 14)
  )

ggsave("sample_map.png", sample_map, width = 8, height = 8)
ggsave("sample_map.pdf", sample_map)
