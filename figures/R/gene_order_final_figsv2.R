library(gridExtra)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(binom)

# Figure 1c
df <- read.table("data/gene_order_final/fig1c.csv", header=TRUE, fill=NA, sep=",")
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

# Figure 2d
df <- read.table("data/gene_order_final/fig2d.csv", header=TRUE, fill=NA, sep=",")
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
df <- read.table("data/gene_order_final/fig3d.csv", header=TRUE, fill=NA, sep=",")
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
df <- read.table("data/gene_order_final/figs2b.csv", header=TRUE, fill=NA, sep=",")
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
df <- read.table("data/gene_order_final/figs4b.csv", header=TRUE, fill=NA, sep=",")
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
df <- read.table("data/gene_order_final/figs4c.csv", header=TRUE, fill=NA, sep=",")
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