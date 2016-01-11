library("ggplot2")
library("dplyr")
library("reshape2")

## Figure 4 - Cases where evolutionarily conserved gene order does not follow 
## assembly order tend to be highly expressed.
df.a <- read.csv("operon_assembly/figures/data/fig4a_ecoli.csv", 
                 header=TRUE, fill=NA)
df.a <- select(df.a, Assembly.same, Assembly.different, Assembly.tied)
colnames(df.a) <- c("Assembly order\nsame as\ngene order",
                    "Assembly order\ndifferent from\ngene order",
                    "Both subunits\nassemble\nsimultaneously")
df.a <- melt(df.a)
df.a <- na.omit(df.a)
fig4a <- ggplot(df.a, aes(variable, value)) +
  geom_boxplot(fill=c("#5BDDFF", "#FD7B80", "#65D760"), outlier.size=0.5) +
  scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                limits=c(0.01, 10000),
                labels=c("0.01", '0.1', '1', '10', '100', '1,000', '10,000')) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))
fig4a
