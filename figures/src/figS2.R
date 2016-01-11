library("ggplot2")
library("dplyr")
library("reshape2")

## Figure S2 - Distribution of intervening genes separating interacting pairs.
incomp1 <-  "operon_assembly/figures/data/figs2a.txt"
therm_comp1 <- "operon_assembly/figures/data/figs2c.txt"
ex_therm_comp1 <- "operon_assembly/figures/data/figs2b.txt"
y2h <- "operon_assembly/figures/data/figs2d.txt"

df.x.ic1 <- read.table(incomp1)
df.x.tc1 <- read.table(therm_comp1)
df.x.etc1 <- read.table(ex_therm_comp1)
df.x.y2h <- read.table(y2h)

# Lots of data reshaping
reshape.figxdata <- function(df){
  colnames(df) <- c("data", "intervening")
  df <- dcast(df, intervening~data)
  df$expected <- df$expected/1000
  df <- melt(df, id.var = "intervening")
  df$variable <- factor(df$variable, 
                        levels = rev(levels(factor(df$variable))))
  df <- rename(df, data = variable)
  return(df)
}

df.x.ic1 <- reshape.figxdata(df.x.ic1)
df.x.tc1 <- reshape.figxdata(df.x.tc1)
df.x.etc1 <- reshape.figxdata(df.x.etc1)
df.x.y2h <- reshape.figxdata(df.x.y2h)

# Generate plots
figx.plot <- function(df){
  figx <- ggplot(df, aes(x = intervening, y=value, fill=data)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Number of intervening genes between interacting pairs") +
    ylab("Number of physically interacting pairs") +
    scale_x_continuous(breaks = seq(0, 14, 2)) +
    theme(text = element_text(size = 6),
          legend.title = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.justification = 'right', 
          legend.position=c(1, 0.9))
  return(figx) 
}

figx.ic1 <- figx.plot(df.x.ic1) +
  scale_y_continuous(breaks = seq(0, 200, 25)) +
  annotate("text", x = 7, y = 190, 
           label = paste(expression(italic("P")), "< 10^-5"),
           size = 4, parse = TRUE)
figx.tc1 <- figx.plot(df.x.tc1) + 
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  annotate("text", x = 4.9, y = 9.5, 
           label = paste(expression(italic("P"))),
           size = 4, parse = TRUE) +
  annotate("text", x = 7.3, y = 9.5, 
           label = "< 0.0004",
           size = 4) # Ugly hack to prevent ggplot converting to e-4
figx.etc1 <- figx.plot(df.x.etc1) + 
  scale_y_continuous(breaks = seq(0, 200, 25)) +
  annotate("text", x = 7, y = 175, 
           label = paste(expression(italic("P")), "< 10^-5"),
           size = 4, parse = TRUE)
figx.y2h <- figx.plot(df.x.y2h) + 
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  annotate("text", x = 7, y = 42, 
           label = paste(expression(italic("P")), "< 0.003"),
           size = 4, parse = TRUE)

# Display plots - save at 3.5 x 3.5 inches
figx.ic1
figx.tc1
figx.etc1
figx.y2h

