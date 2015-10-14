library("gridExtra")
library("ggplot2")
library("dplyr")
library("plyr")

## Figure 1b

data_b <- read.csv("operon_assembly/figures/data/dataset4.csv", header=TRUE, 
                   fill=NA, skip=1)
colnames(data_b) <- c("struc", "sub.A", "sub.B", "gene.A", "gene.B", "species", "interface", 
                      "assembly.order", "gene.fusion", "fusion.conserves", "abundance.A", 
                      "abundance.B", "abundance.A.ecoli", "abundance.B.ecoli", "operon.ID", 
                      "position.A", "position.B", "operon.conservation", "operon.length")
data_b <- filter(data_b, !is.na(abundance.A), !is.na(abundance.B))
combined <- select(data_b, abundance.A, abundance.B)
axis_breaks = c(1, 10, 100, 1000, 10000, 100000)

# Plot panels in figure 1b
fig1b.plot <- function(df){
  ggplot(df, aes(abundance.A, abundance.B)) +
    geom_point() +
    scale_y_log10(breaks=axis_breaks) +
    scale_x_log10(breaks=axis_breaks)
}

fig1b_i <- fig1b.plot(filter(data_b, is.na(operon.ID))) 
fig1b_ii <- fig1b.plot(filter(data_b, !is.na(operon.ID))) 

# Calculate significance of the difference between rho values
rho.difference <- function(df1, df2){
  r1 <- cor.test(df1$abundance.A, df1$abundance.B, method="spearman")
  r2 <- cor.test(df2$abundance.A, df2$abundance.B, method="spearman")
  r1 <- as.numeric(r1[[4]])
  r2 <- as.numeric(r2[[4]])
  return(sqrt((r1-r2)**2))
}

calc.pval <- function(n){
  count <- 0
  real_diff <- rho.difference(filter(data_b, !is.na(operon.ID)),
                              filter(data_b, is.na(operon.ID)))
  for (i in 0:n){
    rand_diff <- rho.difference(sample_n(combined, 89), sample_n(combined, 134))
    if (rand_diff >= real_diff){
      count = count + 1
    }
  }
  return(count/n)
}

## Figure 1c
data_c <- read.csv("operon_assembly/figures/data/fig1c.csv", header=TRUE, fill=NA)
data_c <- select(data_c, op_encoded, diff_tu)
data_c <- melt(data_c)

fig1c <- ggplot(data_c, aes(variable, value)) +
  geom_boxplot(fill=c("darkorange1", "#65D760"), outlier.size=1) +
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))

## Print figures
grid.arrange(fig1b_i, fig1b_ii, ncol=2)
fig1c
