library(ggplot2)
df = data.frame(x=1:10, y=(1:10)+rnorm(10), cal="A")
df2 = data.frame(x=1:10, y=(1:10)*1.5, cal="B")
df = rbind(df, df2)

hom_df = data.frame(hom_value = 5)
inf_cur = data.frame(inf_value = 7)

p = ggplot(df, aes(x=x, y=y, colour=cal)) + geom_line() +
  scale_colour_manual(values = c("A"="red", "B"="blue")) +
  geom_hline(data=hom_df, aes(yintercept=hom_value, linetype="Homogeneous pop."), colour="#008e37") +
  geom_hline(data=inf_cur, aes(yintercept=inf_value, linetype="No switching"), colour="#65009b") +
  scale_linetype_manual(name="Reference models", values=c("Homogeneous pop."="dashed", "No switching"="dotted"),
    guide=guide_legend(override.aes=list(colour=c("#008e37", "#65009b"))))
ggsave("test_legend.pdf", p)
