library(tidyverse)

dstats <- read_tsv("Lo_Cc_Ei_localFstats__10_1.txt")

dstats |>
  ggplot(aes(x = windowStart, y = f_dM)) +
  geom_line() +
  facet_wrap(~chr, scales = "free_x")
