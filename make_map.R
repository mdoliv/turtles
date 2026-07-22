library(rnaturalearth)
library(geobr)
library(tidyverse)
library(ggtext)
library(terra)

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

sample_map <- countries_sa |>
  ggplot() +
  tidyterra::geom_spatraster_rgb(data = ne1_crop) +
  geom_sf(fill = NA) +
  geom_sf(data = states, fill = NA) +
  geom_point(
    data = sample_meta,
    aes(
      x = lon,
      y = lat,
      fill = previous_molecular_id,
      shape = previous_molecular_id
    ),
    size = 4,
    alpha = 0.6
  ) +
  coord_sf(xlim = c(-80, -25), ylim = c(-40, 25), expand = FALSE) +
  #geom_sf_text(data = states, aes(label = postal), size = 3, color = "black") +
  scale_fill_manual(
    name = "Identificação molecular",
    values = c(
      "darkgrey",
      "#F28E2B",
      "#59A14F",
      "#4E79A7",
      "#B07AA1"
    ),
    na.value = "black",
    labels = c(
      "Híbrido",
      "*Ca. caretta*",
      "*Ch. mydas*",
      "*E. imbricata*",
      "*L. olivacea*",
      "Não identificado"
    )
  ) +
  scale_shape_manual(
    name = "Identificação molecular",
    values = c(
      21,
      22,
      23,
      24,
      25
    ),
    na.value = 3,
    labels = c(
      "Híbrido",
      "*Ca. caretta*",
      "*Ch. mydas*",
      "*E. imbricata*",
      "*L. olivacea*",
      "Não identificado"
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    legend.text = element_markdown(size = 12),
    legend.title = element_text(size = 14)
  )

ggsave("sample_map.png", sample_map, width = 8, height = 8)
ggsave("sample_map.pdf", sample_map)
