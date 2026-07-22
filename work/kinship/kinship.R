library(tidyverse)

kinship <- read_tsv("Biallelic.MN10.MX100.MS80.MQ20.kin0")

kinship |>
  ggplot(aes(x = `#IID1`, y = IID2, fill = KINSHIP)) +
  geom_tile() +
  theme_minimal()

kinship |>
  filter(KINSHIP >= 0.25) |>
  write_csv("kinship_first_degree.csv")
