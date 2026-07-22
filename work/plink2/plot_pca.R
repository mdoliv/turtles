library(tidyverse)
library(ggtext)

pca <- read_table(
  "LD_pruned.MN10.MX100.MS80.MQ20.eigenvec",
  col_names = FALSE,
  skip = 1
)
eigenval <- scan("LD_pruned.MN10.MX100.MS80.MQ20.eigenval")

names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

meta <- read_csv("../../data/samples.csv")

pca_with_meta <- pca |>
  left_join(meta, by = c("ind" = "new_id")) |>
  distinct(ind, .keep_all = TRUE)

pve <- data.frame(PC = 1:20, pve = eigenval / sum(eigenval) * 100)

ggplot(pve, aes(as.factor(PC), pve)) +
  geom_col() +
  ylab("Percent variance explained") +
  xlab("PC") +
  theme_minimal()

pca_plot <- ggplot(
  pca_with_meta,
  aes(PC1, PC2, fill = previous_molecular_id, shape = previous_molecular_id)
) +
  geom_point(size = 3) +
  coord_equal() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_minimal() +
  scale_fill_manual(
    name = "Identificação\nmolecular",
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
    name = "Identificação\nmolecular",
    na.value = 3,
    values = c(21, 22, 23, 24, 25),
    labels = c(
      "Híbrido",
      "*Ca. caretta*",
      "*Ch. mydas*",
      "*E. imbricata*",
      "*L. olivacea*",
      "Não identificado"
    )
  ) +
  theme(
    legend.text = element_markdown(),
    #legend.background = element_rect(colour = "black"),
    legend.position = c(0.8, 0.3)
  )
ggsave("pca_plot.pdf", pca_plot)
ggsave("pca_plot.png", pca_plot)

####
pca <- read_table(
  "no_cm.eigenvec",
  col_names = FALSE,
  skip = 1
)
eigenval <- scan("no_cm.eigenval")

names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

meta <- read_csv("../../data/samples.csv")

pca_with_meta <- pca |>
  left_join(meta, by = c("ind" = "new_id")) |>
  distinct(ind, .keep_all = TRUE)

pve <- data.frame(PC = 1:20, pve = eigenval / sum(eigenval) * 100)

ggplot(pve, aes(as.factor(PC), pve)) +
  geom_col() +
  ylab("Percent variance explained") +
  xlab("PC") +
  theme_minimal()

pca_plot_no_cm <- ggplot(
  pca_with_meta,
  aes(PC1, PC2, fill = previous_molecular_id, shape = previous_molecular_id)
) +
  geom_point(size = 3, show.legend = FALSE) +
  coord_equal() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_minimal() +
  scale_fill_manual(
    name = "Identificação molecular",
    values = c(
      "darkgrey",
      "#F28E2B",
      "#4E79A7",
      "#B07AA1"
    ),
    na.value = "black",
    labels = c(
      "Híbrido",
      "*Ca. caretta*",
      "*E. imbricata*",
      "*L. olivacea*",
      "Não identificado"
    )
  ) +
  scale_shape_manual(
    name = "Identificação molecular",
    na.value = 3,
    values = c(21, 23, 24, 25),
    labels = c(
      "Híbrido",
      "*Ca. caretta*",
      "*E. imbricata*",
      "*L. olivacea*",
      "Não identificado"
    )
  ) +
  theme(
    legend.text = element_markdown(),
    legend.background = element_rect(colour = "black"),
    legend.position = c(0.25, 0.15)
  )
ggsave("pca_plot_no_cm.pdf", pca_plot_no_cm)
ggsave("pca_plot_no_cm.png", pca_plot_no_cm)

library(patchwork)

pca_figure <- pca_plot +
  pca_plot_no_cm +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("pca_both.png", pca_figure, width = 8)
