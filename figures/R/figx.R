library("ggplot2")
library("dplyr")
library("reshape2")

## Figure x - Distribution of intervening genes separating interacting pairs.
comp1 <- "operon_assembly/figures/data/intervening_gene_dist_comp1.txt"
excomp1 <-  "operon_assembly/figures/data/intervening_gene_dist_exc_comp1.txt"
incomp1 <-  "operon_assembly/figures/data/intervening_gene_dist_inc_comp1.txt"
y2h <- "operon_assembly/figures/data/intervening_gene_dist_y2h.txt"
df.x.c1 <- read.table(comp1)
df.x.ec1 <- read.table(excomp1)
df.x.ic1 <- read.table(incomp1)
df.x.y2h <- read.table(y2h)

# Lots of data reshaping
reshape.figxdata <- function(df){
  colnames(df) <- c("data", "intervening")
  df <- dcast(df, intervening~data)
  df$expected <- df$expected/1000
  df <- melt(df, id.vars = "intervening")
  df$variable <- factor(df$variable, 
                        levels = rev(levels(factor(df$variable))))
  df <- rename(df, data = variable)
  return(df)
}

df.x.c1 <- reshape.figxdata(df.x.c1)
df.x.ec1 <- reshape.figxdata(df.x.ec1)
df.x.ic1 <- reshape.figxdata(df.x.ic1)
df.x.y2h <- reshape.figxdata(df.x.y2h)

figx.plot <- function(df){
  figx <- ggplot(df, aes(x = intervening, y=value, fill=data)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Number of intervening genes between interacting pairs") +
    ylab("Number of physically interacting pairs") +
    scale_x_continuous(breaks = seq(0, 14, 2)) +
    theme(text = element_text(size = 10),
          legend.title = element_blank())
  return(figx)
}


figx.c1 <- figx.plot(df.x.c1) + 
  scale_y_continuous(breaks = seq(0, 18, 2))
figx.ec1 <- figx.plot(df.x.ec1) + 
  scale_y_continuous(breaks = seq(0, 60, 10))
figx.ic1 <- figx.plot(df.x.ic1) +
  scale_y_continuous(breaks = seq(0, 200, 25))
figx.y2h <- figx.plot(df.x.y2h) + 
  scale_y_continuous(breaks = seq(0, 50, 5))

figx.c1
figx.ic1
figx.ec1
figx.y2h
