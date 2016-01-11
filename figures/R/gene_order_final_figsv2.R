library(gridExtra)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(binom)

## Figure 1b
data1 <- read.csv("operon_assembly/figures/data/dataset1.csv", header=TRUE, 
                  fill=NA, skip=1)
data1 <- rename(data_b,
                abundance.A=Abundance.of.subunit.A.considering.data.from.all.organims, 
                abundance.B=Abundance.of.subunit.B.considering.data.from.all.organisms)
data1 <- filter(data1, !is.na(abundance.A), !is.na(abundance.B))
combined <- select(data1, abundance.A, abundance.B)
axis_breaks = c(1, 10, 100, 1000, 10000, 100000)

# Function to plot panels in figure 1b
fig1b.plot <- function(df){
  ggplot(df, aes(abundance.A, abundance.B)) +
    geom_point() +
    scale_y_log10(breaks=axis_breaks) +
    scale_x_log10(breaks=axis_breaks)
}

fig1b_i <- fig1b.plot(filter(data1, is.na(Operon.ID))) 
fig1b_ii <- fig1b.plot(filter(data1, !is.na(Operon.ID))) 


rho.difference <- function(df1, df2){
  r1 <- cor.test(df1$abundance.A, df1$abundance.B, method="spearman")
  r2 <- cor.test(df2$abundance.A, df2$abundance.B, method="spearman")
  r1 <- as.numeric(r1[[4]])
  r2 <- as.numeric(r2[[4]])
  return(sqrt((r1-r2)**2))
}

calc.pval <- function(n){
  count <- 0
  real_diff <- rho.difference(filter(data1, !is.na(Operon.ID)),
                              filter(data1, is.na(Operon.ID)))
  for (i in 0:n){
    rand_diff <- rho.difference(sample_n(combined, 89), sample_n(combined, 134))
    if (rand_diff >= real_diff){
      count = count + 1
    }
  }
  return(count)
}

# Print figure
grid.arrange(fig1b_i, fig1b_ii, ncol=2)


## Figure 1c
df <- read.csv("operon_assembly/figures/data/fig1c.csv", header=TRUE, fill=NA)
df <- select(df, op_encoded, diff_tu)
df <- rename(df, c("op_encoded"="Operon-encoded\ncomplexes", 
                   "diff_tu"="Complexes\nencoded by different\ntranscriptional units"))
df <- melt(df)
fig1c <- ggplot(df, aes(variable, value)) +
  geom_boxplot(fill=c("darkorange1", "#65D760"), outlier.size=1) +
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))
fig1c

## Figure 2c
df <- read.table('data/prokaryotic_gene_pairs/ecoli_y2h_pairwise.txt', header=TRUE)
df <- mutate(df, int = sqrt((pos1 - pos2)**2) - 1)

build_plt_df <- function(df){
  poss_ppis <- c()
  obs_ppis <- c()
  xvals <- c("0", "1", "2", ">2")
  xvals <- factor(xvals, levels = xvals)
  for (i in 0:3){
    if (i == 3){
      poss_ppis <- append(poss_ppis, length(filter(df, int >= i)$int))
      obs_ppis <- append(obs_ppis, length(filter(df, int >= i, ppi == TRUE)$int))
    } else {
      poss_ppis <- append(poss_ppis, length(filter(df, int == i)$int))
      obs_ppis <- append(obs_ppis, length(filter(df, int == i, ppi == TRUE)$int))      
    }
  }
  plt_df <- data.frame(xvals, poss_ppis, obs_ppis)
  plt_df <- mutate(plt_df, ppi_percs = obs_ppis/poss_ppis*100, 
                   yhi = TRUE, ylo = TRUE, type = TRUE)
  plt_df[plt_df$xvals == "0",]$type <- "Adjacent"
  plt_df[plt_df$xvals != "0",]$type <- "Non-adjacent"
  for (i in 1:4){
    plt_df$ylo[i] <- binom.confint(plt_df$obs_ppis[i], 
                                   plt_df$poss_ppis[i], 0.68)[[5]][5]*100
    plt_df$yhi[i] <- binom.confint(plt_df$obs_ppis[i], 
                                   plt_df$poss_ppis[i], 0.68)[[6]][5]*100
  }
  return(plt_df)
}

plt_data <- function(df){
  plt <- ggplot(df, aes(xvals, ppi_percs, fill=type)) +
    geom_bar(stat="identity", colour = 'black', lwd = 0.4) +
    geom_errorbar(aes(ymax = yhi, ymin=ylo),
                  width=0.1) +
    scale_fill_manual(values = c('firebrick2', 'dodgerblue')) +
    xlab('Intervening genes between\ninteracting pair') +
    ylab('Binary PPIs detected (%)') +
    theme(text = element_text(size = 8),
          legend.key.size = unit(0.25, "cm"),
          legend.justification = 'right', 
          legend.position=c(1,0.85)) +
    annotate("text", x = 1, y = 0.15, label = df$obs_ppis[1], size = 2.5) +
    annotate("text", x = 2, y = 0.15, label = df$obs_ppis[2], size = 2.5) +
    annotate("text", x = 3, y = 0.15, label = df$obs_ppis[3], size = 2.5) +
    annotate("text", x = 4, y = 0.15, label = df$obs_ppis[4], size = 2.5) +
    guides(fill=guide_legend(title="Gene pairs"))
  plt
}
plt_df <- build_plt_df(df)
plt_data(plt_df)

altPlotter <- function(df){
  xlevs <- c("0", "1", "2", "3", "4", "5", "6", ">6")
  #   df <- filter(df, ppi == TRUE)
  df <- mutate(df, type = "Non-adjacent")
  df[df$int == 0,]$type <- "Adjacent"
  #   df[df$int > 6,]$int <- ">6"
  df$int <- factor(df$int, levels = xlevs)
  plt <- ggplot(df) +
    geom_bar(aes(int, fill = ppi))
  return(plt)
}
altPlotter(df)

# Figure 2d
df <- read.table("operon_assembly/figures/data/fig2d.csv", header=TRUE, fill=NA, sep=",")
df <- rename(df, c("Adjacent"="Adjacent\ngenes",
                   "Non.adjacent"="Non-adjacent\ngenes", 
                   "DiffTU"="Subunits\nencoded by different\ntranscription units"))
df <- melt(df)
fig2d <- ggplot(df, aes(variable, value)) +
  geom_boxplot(fill=c("darkorange1", "skyblue3", "#65D760"), outlier.size=1) +
  scale_y_log10(breaks=c(200, 500, 1000, 2000, 5000, 10000, 20000)) +
  ylab(expression(paste("Interface size (\uc5"^"2",")"))) +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))
fig2d

# Figure 3d
df <- read.table("operon_assembly/figures/data/fig3d.csv", header=TRUE, fill=NA, sep=",")
df <- select(df, Assembly.same, Assembly.different, Assembly.tied)
df <- rename(df, c("Assembly.same"="Assembly order\nsame as\ngene order", 
                   "Assembly.different"="Assembly order\ndifferent from\ngene order", 
                   "Assembly.tied"="Both subunits\nassemble\nsimultaneously"))
df <- melt(df)
fig3d <- ggplot(df, aes(variable, value)) +
  geom_boxplot(fill=c("#5BDDFF", "#FD7B80", "#65D760"), outlier.size=1) +
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))
fig3d

# Figure s2b
df <- read.table("operon_assembly/figures/data/figs2b.csv", header=TRUE, fill=NA, sep=",")
df <- select(df, Assembly.same, Assembly.different, Assembly.tied)
df <- rename(df, c("Assembly.same"="Assembly order\nsame as\ngene order", 
                   "Assembly.different"="Assembly order\ndifferent from\ngene order", 
                   "Assembly.tied"="Both subunits\nassemble\nsimultaneously"))
df <- melt(df)
figs2b <- ggplot(df, aes(variable, value)) +
  geom_boxplot(fill=c("#5BDDFF", "#FD7B80", "#65D760"), outlier.size=1) +
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))
figs2b

# Figure s4b
df <- read.table("operon_assembly/figures/data/figs4b.csv", header=TRUE, fill=NA, sep=",")
df <- select(df, Op.encoded, Diff.TU)
df <- rename(df, c("Op.encoded"="Operon-encoded\ncomplexes", 
                   "Diff.TU"="Complexes\nencoded by different\ntranscriptional units"))
df <- melt(df)
figs4b <- ggplot(df, aes(variable, value)) +
  geom_boxplot(fill=c("darkorange1", "#65D760"), outlier.size=1) +
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))
figs4b

# Figure s4c
df <- read.table("operon_assembly/figures/data/figs4c.csv", header=TRUE, fill=NA, sep=",")
df <- select(df, Assembly.same, Assembly.different, Assembly.tied)
df <- rename(df, c("Assembly.same"="Assembly order\nsame as\ngene order", 
                   "Assembly.different"="Assembly order\ndifferent from\ngene order", 
                   "Assembly.tied"="Both subunits\nassemble\nsimultaneously"))
df <- melt(df)
figs4c <- ggplot(df, aes(variable, value)) +
  geom_boxplot(fill=c("#5BDDFF", "#FD7B80", "#65D760"), outlier.size=1) +
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) +
  ylab("Protein abundance (ppm)") +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color='black'))
figs4c