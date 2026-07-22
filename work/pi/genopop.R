library(GenoPop)
library(tidyverse)

# Total number of sites in the VCF file
total_seqlen <- 1065048

lo_pi <- Pi("Lo.no_related.vcf.gz", total_seqlen, threads = 16)
cc_br_pi <- Pi("Cc_BR.All.vcf.gz", total_seqlen, threads = 16)
cc_erg_pi <- Pi("Cc_ERG.All.vcf.gz", total_seqlen, threads = 16)
cm_pi <- Pi("Cm.All.vcf.gz", total_seqlen, threads = 16)
ei_pi <- Pi("Ei.All.vcf.gz", total_seqlen, threads = 16)

lo_het <- Heterozygosity("Lo.no_related.vcf.gz", threads = 16)
cc_br_het <- Heterozygosity("Cc_BR.All.vcf.gz", threads = 16)
cc_erg_het <- Heterozygosity("Cc_ERG.All.vcf.gz", threads = 16)
cm_het <- Heterozygosity("Cm.All.vcf.gz", threads = 16)
ei_het <- Heterozygosity("Ei.All.vcf.gz", threads = 16)

lo_theta <- WattersonsTheta("Lo.no_related.vcf.gz", total_seqlen, threads = 16)
cc_br_theta <- WattersonsTheta("Cc_BR.All.vcf.gz", total_seqlen, threads = 16)
cc_erg_theta <- WattersonsTheta("Cc_ERG.All.vcf.gz", total_seqlen, threads = 16)
cm_theta <- WattersonsTheta("Cm.All.vcf.gz", total_seqlen, threads = 16)
ei_theta <- WattersonsTheta("Ei.All.vcf.gz", total_seqlen, threads = 16)

lo_tajima <- TajimasD("Lo.no_related.vcf.gz", total_seqlen, threads = 16)
cc_br_tajima <- TajimasD("Cc_BR.All.vcf.gz", total_seqlen, threads = 16)
cc_erg_tajima <- TajimasD("Cc_ERG.All.vcf.gz", total_seqlen, threads = 16)
cm_tajima <- TajimasD("Cm.All.vcf.gz", total_seqlen, threads = 16)
ei_tajima <- TajimasD("Ei.All.vcf.gz", total_seqlen, threads = 16)

summary_stats <- data.frame(
  pop = rep(c("Lo", "Cc BR", "Cc ERG", "Cm", "Ei"), 4),
  stat = c(
    rep("Pi", 5),
    rep("Heterozygosity", 5),
    rep("Theta", 5),
    rep("Tajima's D", 5)
  ),
  value = c(
    lo_pi,
    cc_br_pi,
    cc_erg_pi,
    cm_pi,
    ei_pi,
    lo_het,
    cc_br_het,
    cc_erg_het,
    cm_het,
    ei_het,
    lo_theta,
    cc_br_theta,
    cc_erg_theta,
    cm_theta,
    ei_theta,
    lo_tajima,
    cc_br_tajima,
    cc_erg_tajima,
    cm_tajima,
    ei_tajima
  )
)
write_csv(summary_stats, "genopop_summary.csv")

summary_stats_wide <- summary_stats |>
  pivot_wider(names_from = stat, values_from = value)
write_csv(summary_stats_wide, "genopop_summary_wide.csv")

summary_stats <- read_csv("genopop_summary.csv")

diversity_plot <- summary_stats |>
  filter(stat != "Tajima's D") |>
  mutate(pop = fct_reorder(pop, value)) |>
  ggplot(aes(y = pop, x = value, fill = pop)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(
    ~stat,
    nrow = 3
  ) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(expand = NA) +
  scale_fill_manual(
    values = c(
      "#377eb8",
      "#ff7f00",
      "#984ea3",
      "#e41a1c",
      "#4daf4a"
    )
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggsave("diversity_stats.png", diversity_plot)

# Fst calculations
cc_br_samples <- read_lines("Cc_BR_samples.txt")
cc_erg_samples <- read_lines("Cc_ERG_samples.txt")
lo_samples <- read_lines("Lo_samples_no_related.txt")
ei_samples <- read_lines("Ei_samples.txt")
cm_samples <- read_lines("Cm_samples.txt")

cc_all_fst <- Fst(
  "Cc_All.All.vcf.gz",
  cc_br_samples,
  cc_erg_samples,
  threads = 16
)
cc_ei_fst <- Fst(
  "../dadi/All.cc_ei.vcf.gz",
  cc_br_samples,
  ei_samples,
  threads = 16
)
cc_lo_fst <- Fst(
  "../dadi/All.cc_lo.vcf.gz",
  cc_br_samples,
  lo_samples,
  threads = 16
)
cc_cm_fst <- Fst(
  "../dadi/All.cc_cm.vcf.gz",
  cc_br_samples,
  cm_samples,
  threads = 16
)
ei_lo_fst <- Fst(
  "../dadi/All.ei_lo.vcf.gz",
  ei_samples,
  lo_samples,
  threads = 16
)
ei_cm_fst <- Fst(
  "../dadi/All.ei_cm.vcf.gz",
  ei_samples,
  cm_samples,
  threads = 16
)
cm_lo_fst <- Fst(
  "../dadi/All.cm_lo.vcf.gz",
  cm_samples,
  lo_samples,
  threads = 16
)

fst_data <- data.frame(
  pop_1 = c(rep("Cc BR", 4), "Ei", "Ei", "Cm"),
  pop_2 = c("Cc ERG", "Ei", "Lo", "Cm", "Lo", "Cm", "Lo"),
  value = c(
    cc_all_fst,
    cc_ei_fst,
    cc_lo_fst,
    cc_cm_fst,
    ei_lo_fst,
    ei_cm_fst,
    cm_lo_fst
  )
)
write_csv(fst_data, "genopop_fst.csv")

fst_data |>
  ggplot(aes(x = pop_1, y = pop_2, fill = value)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_viridis_c()

# Dxy calculations
cc_all_dxy <- Dxy(
  "Cc_All.All.vcf.gz",
  cc_br_samples,
  cc_erg_samples,
  total_seqlen,
  threads = 16
)
cc_ei_dxy <- Dxy(
  "../dadi/All.cc_ei.vcf.gz",
  cc_br_samples,
  ei_samples,
  total_seqlen,
  threads = 16
)
cc_lo_dxy <- Dxy(
  "../dadi/All.cc_lo.vcf.gz",
  cc_br_samples,
  lo_samples,
  total_seqlen,
  threads = 16
)
cc_cm_dxy <- Dxy(
  "../dadi/All.cc_cm.vcf.gz",
  cc_br_samples,
  cm_samples,
  total_seqlen,
  threads = 16
)
ei_lo_dxy <- Dxy(
  "../dadi/All.ei_lo.vcf.gz",
  ei_samples,
  lo_samples,
  total_seqlen,
  threads = 16
)
ei_cm_dxy <- Dxy(
  "../dadi/All.ei_cm.vcf.gz",
  ei_samples,
  cm_samples,
  total_seqlen,
  threads = 16
)
cm_lo_dxy <- Dxy(
  "../dadi/All.cm_lo.vcf.gz",
  cm_samples,
  lo_samples,
  total_seqlen,
  threads = 16
)

dxy_data <- data.frame(
  pop_1 = c(rep("Cc BR", 4), "Ei", "Ei", "Cm"),
  pop_2 = c("Cc ERG", "Ei", "Lo", "Cm", "Lo", "Cm", "Lo"),
  value = c(
    cc_all_dxy,
    cc_ei_dxy,
    cc_lo_dxy,
    cc_cm_dxy,
    ei_lo_dxy,
    ei_cm_dxy,
    cm_lo_dxy
  )
)
write_csv(dxy_data, "genopop_dxy.csv")
