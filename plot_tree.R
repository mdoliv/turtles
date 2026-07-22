library(treedataverse)
library(stringr)

turtles_tree <- read.newick("data/turtle_species_renamed.nwk")
turtles_tree$tip.label <- str_replace(turtles_tree$tip.label, "_", " ")

tree_plot <- ggtree(turtles_tree, branch.length = "none") +
  geom_tiplab(fontface = 3) +
  geom_rootedge(0.25) +
  scale_x_continuous(expand = c(0, 3))

ggsave("tree_plot.png", tree_plot)
ggsave("tree_plot.pdf", tree_plot)
