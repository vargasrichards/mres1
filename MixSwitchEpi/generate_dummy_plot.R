library(ggplot2)
p <- ggplot() +
  geom_hline(aes(yintercept=1, linetype="Homogeneous pop."), color="#E41A1C", size=1) +
  geom_hline(aes(yintercept=2, linetype="Heterogeneous pop. \n no switching"), color="#4DAF4A", size=1) +
  scale_linetype_manual(
    name = "Reference models",
    values = c("Homogeneous pop." = "dashed", "Heterogeneous pop. \n no switching" = "dotted"),
    guide = guide_legend(override.aes = list(colour = c("#4DAF4A", "#E41A1C")))
  )
ggsave("dummy_output.png", p)
