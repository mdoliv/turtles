library(tidyverse)

ddrad <- read_tsv(
  "ddRAD.insert.tsv",
  col_names = c("size", "total", "inward", "outward", "other")
) |>
  mutate(library = "ddRAD")
rad3 <- read_tsv(
  "3RAD.insert.tsv",
  col_names = c("size", "total", "inward", "outward", "other")
) |>
  mutate(library = "3RAD")

all_rad <- rbind(ddrad, rad3)

insert_plot <- all_rad |>
  ggplot(aes(x = size, y = inward, fill = library, colour = library)) +
  geom_area(alpha = 0.5) +
  scale_fill_discrete(name = "Library") +
  scale_colour_discrete(name = "Library") +
  labs(x = "Insert size (bp)", y = "Inward oriented pairs") +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

ggsave("insert_size.png", insert_plot)
