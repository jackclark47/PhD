library(openxlsx)
library(ggplot2)

data <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/significant_loci_matches.xlsx")
data = out_main
#Plot NoVar
ggplot(data, aes(x = as.numeric(exp), y = percent_both, color = sig_change)) +
  geom_point() +
  xlab("log2 fold change (MAX/MIN)") +
  ylab("Percentage match") +
  guides(fill=guide_legend(title="New Legend Title")) +
  theme(legend.title=element_blank(),
        axis.text = element_text(size=20),
        axis.title=element_text(size=25),
        legend.text = element_text(size=20))
#Plot N188Var
ggplot(data, aes(x = as.numeric(exp), y = percent_both, color = sig_change)) +
  geom_point() +
  xlab("log2 fold change (MAX/MIN)") +
  ylab("Percentage match") +
  guides(color=guide_legend(title="Significant change"))

#Plot Multiple Var
ggplot(data, aes(x = as.numeric(exp), y = percent_both, color = sig_change)) +
  geom_point() +
  xlab("log2 fold change (MAX/MIN)") +
  ylab("Percentage match") +
  guides(color=guide_legend(title="New Legend Title"))
