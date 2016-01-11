library("grid")
library("reshape2")
library("ggplot2")
library("plyr")
library("dplyr")
library("binom")

## Figure 2b
df.b <- read.csv("operon_assembly/figures/data/fig2b.csv", header=TRUE, 
                 fill=NA, skip=1)
colnames(df.b) <- c("struc", "sub.A", "sub.B", "gene.A", "gene.B", "species", 
                    "interface", "assembly.order", "gene.fusion", 
                    "fusion.conserves", "abundance.A", "abundance.B", 
                    "abundance.A.ecoli", "abundance.B.ecoli", "operon.ID", 
                    "position.A", "position.B", "operon.conservation", 
                    "operon.length")
df.b <- filter(df.b, !is.na(operon.ID))
df.b <- mutate(df.b, intervening = sqrt((position.A - position.B)**2) - 1, 
                 type = "Non-Adjacent", int_existence = "Yes",
                 col = "darkorange")
df.b[df.b$intervening == 0,]$type <- "Adjacent"
df.b[df.b$interface <= 200,]$int_existence <- "No"
df.b[df.b$int_existence == "No",]$col <- "darkorange2"
df.b[df.b$intervening > 6,]$intervening <- ">6"
df.b$type <- factor(df.b$type, levels = c("Adjacent", "Non-adjacent"))
df.b$int_existence <- factor(df.b$int_existence, levels = c("Yes", "No"))

# Returns stacked barplot
stacked.plot <- function(df){
  xlevs <- c("0", "1", "2", "3", "4", "5", "6", "", ">6")
  df$intervening <- factor(df$intervening, levels = xlevs)
  plt <- ggplot() +
    geom_bar(data = filter(df, intervening == 0), 
             aes(intervening, fill = col)) +
    geom_bar(data = filter(df, intervening != 0), 
             aes(intervening, fill = int_existence)) +
    scale_fill_manual(values = c("darkorange1", "goldenrod1", 
                                 "skyblue1", "deepskyblue3")) +
    scale_x_discrete(limits = xlevs) +
    theme(legend.position = "none",
          text = element_text(size = 8)) +
    xlab("Intervening genes between pair") +
    ylab("Gene pairs") +
    geom_vline(xintercept = c(8), linetype = 2, alpha = 0.5, lwd = 0.6)
  return(plt)
}

# Create plots
fig2b <- stacked.plot(df.b)
fig2b_complex1 <- stacked.plot(filter(df.b, operon.ID == 116617 | 
                                          operon.ID == 580994))
fig2b_not_complex1 <- stacked.plot(filter(df.b, operon.ID != 116617 &
                                              operon.ID != 580994 ))
fig2b_3to10 <- stacked.plot(filter(df.b, operon.length > 2, 
                                     operon.length <= 10))
fig2b_10plus <- stacked.plot(filter(df.b, operon.length > 9))

## Figure 2c
df.c <- read.table("operon_assembly/figures/data/fig2c.txt", 
                   header=TRUE)
df.c <- mutate(df.c, int = sqrt((position_A - position_B)**2) - 1)

# Format data for barplot
build.plotdata <- function(df){
  poss_ppis <- c()
  obs_ppis <- c()
  xvals <- c("0", "1", "2", ">2")
  xvals <- factor(xvals, levels = xvals)
  for (i in 0:3){
    if (i == 3){
      poss_ppis <- append(poss_ppis, length(filter(df, int >= i)$int))
      obs_ppis <- append(obs_ppis, 
                         length(filter(df, int >= i,
                                       y2h_ppi_detected == TRUE)$int))
    } else {
      poss_ppis <- append(poss_ppis, length(filter(df, int == i)$int))
      obs_ppis <- append(obs_ppis, 
                         length(filter(df, int == i,
                                       y2h_ppi_detected == TRUE)$int))      
    }
  }
  plt_df <- data.frame(xvals, poss_ppis, obs_ppis)
  plt_df <- mutate(plt_df, ppi_percs = obs_ppis/poss_ppis*100, 
                   yhi = TRUE, ylo = TRUE, type = TRUE)
  plt_df[plt_df$xvals == "0",]$type <- "Adjacent"
  plt_df[plt_df$xvals != "0",]$type <- "Non-adjacent"
  for (i in 1:4){
    plt_df$ylo[i] <- binom.confint(plt_df$obs_ppis[i], 
                                   plt_df$poss_ppis[i], 0.68)[[5]][11]*100
    plt_df$yhi[i] <- binom.confint(plt_df$obs_ppis[i], 
                                   plt_df$poss_ppis[i], 0.68)[[6]][11]*100
  }
  return(plt_df)
}

# Builds plot
df.c <- build.plotdata(df.c)
fig2c <- ggplot(df.c, aes(xvals, ppi_percs, fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax = yhi, ymin=ylo),
                width=0.1, alpha=0.75) +
  scale_fill_manual(values = c("darkorange", "skyblue3")) +
  xlab("Intervening genes") +
  ylab("% of E.coli gene pairs with binary\n protein-protein interactions detected") +
  theme(text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.justification = "right", 
        legend.position=c(1,0.85)) +
  annotate("text", x = 1, y = 0.15, 
           label = paste(df.c$obs_ppis[1], df.c$poss_ppis[1], sep="/"), 
           size = 3, alpha=0.75) +
  annotate("text", x = 2, y = 0.15, 
           label = paste(df.c$obs_ppis[2], df.c$poss_ppis[2], sep="/"), 
           size = 3, alpha=0.75) +
  annotate("text", x = 3, y = 0.15, 
           label = paste(df.c$obs_ppis[3], df.c$poss_ppis[3], sep="/"), 
           size = 3, alpha=0.75) +
  annotate("text", x = 4, y = 0.15, 
           label = paste(df.c$obs_ppis[4], df.c$poss_ppis[4], sep="/"), 
           size = 3, alpha=0.75) +
  guides(fill=guide_legend(title="Gene pairs"))
fig2c

# Figure 2d
df.d <- read.csv("operon_assembly/figures/data/fig2d.csv", header=TRUE, fill=NA)
colnames(df.d) <- c("Adjacent\ngenes", "Non-adjacent\ngenes",
                    "Subunits encoded\nby different\ntranscriptional units")
df.d <- melt(df.d)

fig2d <- ggplot(df.d, aes(variable, value)) +
  geom_boxplot(fill=c("darkorange1", "skyblue3", "#65D760"), outlier.size=0.5) +
  scale_y_log10(breaks=c(200, 500, 1000, 2000, 5000, 10000, 20000)) +
  ylab(expression(paste("Interface size (\uc5"^"2",")"))) +
  theme(text = element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(color="black"))

## Display plots
fig2b
fig2c
fig2d
