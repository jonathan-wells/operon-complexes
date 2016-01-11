library("gridExtra")
library("ggplot2")
library("dplyr")
library("reshape2")
library("scales")

## Figure 1b
# To plot figure S1a, use abundance.A/B instead of abundance.A/B.ecoli
# For figure S1b use dataset1_weiss and abundance.A/B. Will probably need to
# change scales as well.
df.b <- read.csv("operon_assembly/figures/data/dataset1_newpaxdb.csv", 
                 header=TRUE, fill=NA, skip=1)
colnames(df.b) <- c("struc", "sub.A", "sub.B", "gene.A", "gene.B", "species", 
                    "interface", "assembly.order", "gene.fusion", 
                    "fusion.conserves", "abundance.A", "abundance.B", 
                    "abundance.A.ecoli", "abundance.B.ecoli", "operon.ID", 
                    "position.A", "position.B", "operon.conservation")
df.b <- filter(df.b, !is.na(abundance.A.ecoli), !is.na(abundance.B.ecoli))

# Plot panels in figure 1b
fig1b.plot <- function(df, pcolour){
  axis_breaks <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
  axis_labels <- c("0.01", "0.1", "1", "10", "100", "1,000", "10,000")
  ggplot(df, aes(abundance.A.ecoli, abundance.B.ecoli)) +
    geom_point(colour=pcolour) +
    geom_point(pch=21) +
    ylab("First protein abundance (ppm)") +
    xlab("Second protein abundance (ppm)") +
    scale_y_log10(breaks=axis_breaks, limits=c(0.01, 20000),
                  labels=axis_labels) +
    scale_x_log10(breaks=axis_breaks, limits=c(0.01, 50000),
                  labels=axis_labels) +
    theme(text = element_text(size=10),
          axis.text=element_text(color='black'))
}

fig1b_i <- fig1b.plot(filter(df.b, is.na(operon.ID)), '#65D760') 
fig1b_ii <- fig1b.plot(filter(df.b, !is.na(operon.ID)), 'darkorange1')

# Calculate significance of the difference between rho values
rho.difference <- function(df1, df2){
  r1 <- cor.test(df1$abundance.A.ecoli, df1$abundance.B.ecoli, 
                 method="spearman")
  r2 <- cor.test(df2$abundance.A.ecoli, df2$abundance.B.ecoli, 
                 method="spearman")
  r1 <- as.numeric(r1[[4]])
  r2 <- as.numeric(r2[[4]])
  return(c(sqrt((r1-r2)**2), r1, r2))
}

calc.pval <- function(n){
  combined = select(df.b, abundance.A.ecoli, abundance.B.ecoli)
  count <- 0
  real_diff <- rho.difference(filter(df.b, !is.na(operon.ID)),
                              filter(df.b, is.na(operon.ID)))
  for (i in 0:n){
    rand_diff <- rho.difference(sample_n(combined, 89), 
                                sample_n(combined, 134))
    if (rand_diff[1] >= real_diff[1] & rand_diff[2] >= rand_diff[3]){
      count = count + 1
    }
  }
  return(count/n)
}


## Figure 1c
# For figures S1a and b use fig1c_all.csv or fig1c_weismann.csv and adjust scale
df.c <- read.csv("operon_assembly/figures/data/fig1c_ecoli.csv", 
                 header=TRUE, fill=NA)
df.c <- select(df.c, op_encoded, diff_tu)
col <- c("Complexes\nencoded by different\ntranscriptional units",
         "Operon-encoded\ncomplexes")
colnames(df.c) <- col
df.c <- melt(df.c)
df.c$variable <- factor(df.c$variable, levels = levels(rev(factor(col))))
df.c <- na.omit(df.c)
fig1c <- ggplot(df.c, aes(variable, value)) +
  geom_boxplot(fill=c("#65D760", "darkorange1"), outlier.size=0.5) +
  scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100, 1000, 10000), 
                limits=c(0.01, 10000), 
                labels=c("0.01", "0.1", "1" ,"10", "100", 
                         "1,000", "10,000")) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))

## Display plots
grid.arrange(fig1b_i, fig1b_ii, ncol=2)
fig1c
