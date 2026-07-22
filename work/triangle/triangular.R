library(triangulaR)
library(vcfR)
library(tidyverse)
library(patchwork)
library(ggtext)

cc_ei_vcf <- read.vcfR("cc_ei_hybrids.vcf.gz")
cc_ei_popmap <- read_tsv("cc_ei_hybrids.popmap.tsv", col_names = c("id", "pop"))

cc_ei_diff <- alleleFreqDiff(
  cc_ei_vcf,
  cc_ei_popmap,
  p1 = "Cc",
  p2 = "Ei",
  difference = 1
)
cc_ei_hi <- hybridIndex(cc_ei_diff, cc_ei_popmap, p1 = "Cc", p2 = "Ei")
cc_ei_plot <- triangle.plot(cc_ei_hi, cex = 4) +
  scale_color_manual(
    name = "ADMIXTURE assignment",
    values = c(
      "#F28E2B",
      "#E15759",
      "grey56"
    ),
    labels = c("*Ca. caretta*", "*E. imbricata*", "Hybrid")
  ) +
  labs(x = "Hybrid index", y = "Interclass heterozygosity") +
  theme(
    legend.background = element_rect(colour = "black"),
    legend.text = element_markdown()
  )

ggsave("cc_ei_hybrids.png", cc_ei_plot, width = 8, height = 6.67)

#cc_ei_missing <- missing.plot(cc_ei_hi, cex = 4) +
#  scale_colour_viridis_c(name = "Fraction missing") +
#  labs(title = "Cc x Ei")
#
#ggsave("cc_ei_missing.png", cc_ei_missing, width = 8, height = 6.67)
#ggsave("cc_ei_missing.pdf", cc_ei_missing, width = 8, height = 6.67)

###

cc_lo_vcf <- read.vcfR("cc_lo_hybrids.vcf.gz")
cc_lo_popmap <- read_tsv("cc_lo_hybrids.popmap.tsv", col_names = c("id", "pop"))

cc_lo_diff <- alleleFreqDiff(
  cc_lo_vcf,
  cc_lo_popmap,
  p1 = "Cc",
  p2 = "Lo",
  difference = 1
)
cc_lo_hi <- hybridIndex(cc_lo_diff, cc_lo_popmap, p1 = "Cc", p2 = "Lo")
cc_lo_plot <- triangle.plot(cc_lo_hi, cex = 4) +
  scale_color_manual(
    name = "ADMIXTURE assignment",
    values = c(
      "#F28E2B",
      "grey56",
      "#B07AA1"
    ),
    labels = c("*Ca. caretta*", "Hybrid", "*L. olivacea*")
  ) +
  labs(x = "Hybrid index", y = "Interclass heterozygosity") +
  theme(
    legend.background = element_rect(colour = "black"),
    legend.text = element_markdown()
  )

ggsave("cc_lo_hybrids.png", cc_lo_plot, width = 8, height = 6.67)

#cc_lo_missing <- missing.plot(cc_lo_hi, cex = 4) +
#  scale_colour_viridis_c(name = "Fraction missing") +
#  labs(title = "Cc x Lo")
#
#ggsave("cc_lo_missing.png", cc_lo_missing, width = 8, height = 6.67)
#ggsave("cc_lo_missing.pdf", cc_lo_missing, width = 8, height = 6.67)

###

ei_lo_vcf <- read.vcfR("ei_lo_hybrids.vcf.gz")
ei_lo_popmap <- read_tsv("ei_lo_hybrids.popmap.tsv", col_names = c("id", "pop"))

ei_lo_diff <- alleleFreqDiff(
  ei_lo_vcf,
  ei_lo_popmap,
  p1 = "Ei",
  p2 = "Lo",
  difference = 1
)
ei_lo_hi <- hybridIndex(ei_lo_diff, ei_lo_popmap, p1 = "Ei", p2 = "Lo")
ei_lo_plot <- triangle.plot(ei_lo_hi, cex = 4) +
  scale_color_manual(
    name = "ADMIXTURE assignment",
    values = c(
      "#E15759",
      "grey56",
      "#B07AA1"
    ),
    labels = c("*E. imbricata*", "Hybrid", "*L. olivacea*")
  ) +
  labs(x = "Hybrid index", y = "Interclass heterozygosity") +
  theme(
    legend.background = element_rect(colour = "black"),
    legend.text = element_markdown()
  )

ggsave("ei_lo_hybrids.png", ei_lo_plot, width = 8, height = 6.67)

#ei_lo_missing <- missing.plot(ei_lo_hi, cex = 4) +
#  scale_colour_viridis_c(name = "Fraction missing") +
#  labs(title = "Ei x Lo")
#
#ggsave("ei_lo_missing.png", ei_lo_missing, width = 8, height = 6.67)
#ggsave("ei_lo_missing.pdf", ei_lo_missing, width = 8, height = 6.67)

###

composite_plot <- (cc_ei_plot + cc_lo_plot) /
  (ei_lo_plot + plot_spacer()) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("composite_plot.pdf", composite_plot, width = 8, height = 6.67)

###
#composite_missing <- (cc_ei_missing + cc_lo_missing) /
#  (ei_lo_missing + plot_spacer())
#
#ggsave("composite_missing.png", composite_missing, width = 10, height = 6.67)
###
write_csv(cc_ei_hi, "cc_ei_hi.csv")
write_csv(cc_lo_hi, "cc_lo_hi.csv")
write_csv(ei_lo_hi, "ei_lo_hi.csv")
