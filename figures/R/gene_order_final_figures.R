library(gridExtra)
library(reshape2)
library(ggplot2)
library(dplyr)
library(binom)

#############################################################################
# Abundance vs. Assembly order matches & same operon diff tu
pa_all <- read.table('/Users/jonwells/Downloads/panel_a_all.csv', 
                     header=TRUE, fill=NA)
# pb_all <- read.table('/Users/jonwells/Downloads/panel_b_all.csv', 
#                      header=TRUE, fill=NA, sep=',')
pa_ecoli <- read.table('/Users/jonwells/Downloads/panel_a_ecoli.csv', 
                     header=TRUE, fill=NA)
# pb_ecoli <- read.table('/Users/jonwells/Downloads/panel_b_ecoli.csv', 
#                      header=TRUE, fill=NA, sep=',')
colnames(pa_all)[1] <- 'Assembly order same\nas gene order'
colnames(pa_all)[2] <- 'Assembly order different\nfrom gene order'
# colnames(pb_all)[1] <- 'Same-operon subunits\n'
# colnames(pb_all)[2] <- 'Different TU subunits\n'
colnames(pa_ecoli)[1] <- 'Assembly order same\nas gene order'
colnames(pa_ecoli)[2] <- 'Assembly order different\nfrom gene order'
# colnames(pb_ecoli)[1] <- 'Same-operon subunits\n'
# colnames(pb_ecoli)[2] <- 'Different TU subunits\n'

pa_all <- melt(pa_all)
# pb_all <- melt(pb_all)
pa_ecoli <- melt(pa_ecoli)
# pb_ecoli <- melt(pb_ecoli)

pa.plotter <- function(df, y_ax=TRUE){
  plt <- ggplot(df, aes(variable, value)) +
    geom_boxplot(fill = c('cyan3', 'darkorchid'), outlier.size = 1)  + 
    scale_y_log10(breaks=c(1,10,100,1000,10000), limits=c(0.5, 30000)) +
#     scale_y_log10() +
    ylab("Abundance (ppm)") +
    theme(text = element_text(size=8),
          axis.title.x=element_blank())
  if (y_ax == FALSE){
    plt <- plt + theme(axis.title.y=element_blank())
  }
  return(plt)
}

pb.plotter <- function(df, y_ax=TRUE){
  plt <- ggplot(df, aes(variable, value)) +
    geom_boxplot(fill = c('cadetblue3', 'darkolivegreen3'), outlier.size = 1) + 
    scale_y_log10(breaks=c(1,10,100,1000,10000), limits=c(0.18, 50000)) +
    ylab("Abundance (ppm)") +
    theme(text = element_text(size=8),
          axis.title.x=element_blank())
  if (y_ax == FALSE){
    plt <- plt + theme(axis.title.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.text.y=element_blank())
  }
  return(plt)
}

pa_all_plot <- pa.plotter(pa_all)
pa_ecoli_plot <- pa.plotter(pa_ecoli)
# pb_all_plot <- pb.plotter(pb_all)
# pb_ecoli_plot <- pb.plotter(pb_ecoli)

# grid.arrange(pa_all_plot, pb_all_plot, ncol=2)
# grid.arrange(pa_ecoli_plot, pb_ecoli_plot, ncol=2)
grid.arrange(pa_all_plot, pa_ecoli_plot, ncol=2)

#############################################################################
# Interface type and intervening genes vs interface

gene_pairs <- read.table('data/prokaryotic_gene_pairs/pairs_w_dimers_corrected.txt',
                         header=TRUE, fill=NA)
gene_pairs <- mutate(gene_pairs, intervening = TRUE)
gene_pairs <- mutate(gene_pairs, type = TRUE)
gene_pairs[is.na(gene_pairs$opID),]$type <- "Different TU"
gene_pairs[is.na(gene_pairs$opID),]$intervening <- NA
gene_pairs$intervening <- sqrt((gene_pairs$op2 - gene_pairs$op1)**2)-1
gene_pairs[!is.na(gene_pairs$opID) & 
             gene_pairs$intervening == 0,]$type <- "Adjacent"
gene_pairs[!is.na(gene_pairs$opID) & 
             gene_pairs$intervening != 0,]$type <- "Non-adjacent"
ops_only <- na.omit(gene_pairs)
comp1 <- filter(gene_pairs, opID == 116617)
gene_pairs$type <- factor(gene_pairs$type, 
                          levels = c("Adjacent", "Non-adjacent", 'Different TU'))

type_vs_interface <- function(df){
  # Returns boxplot of interface size vs pos. of gene in or out of operon.
  plt <- ggplot(df[df$int >= 200,], aes(type, int)) + 
    geom_boxplot(fill = c('darkorange', 'skyblue3', '#C22326'), outlier.size = 1) +
    xlab('Gene pairs\n') + ylab(expression(paste("Interface size (\uc5"^"2",")"))) +
    scale_y_log10(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000)) +
    theme(text = element_text(size=8), 
          axis.title.x = element_text(vjust = -0.5)) +
#               axis.title.y = element_text(vjust = 0.5)) +
    annotate("text", x = 1.2, y = 500, 
             label = length(df[df$type == "Adjacent" & df$int >= 200,]$int), size = 2.5) +
    annotate("text", x = 2.2, y = 500, 
             label = length(df[df$type == "Non-adjacent" & df$int >= 200,]$int), size = 2.5) +
    annotate("text", x = 3.2, y = 500, 
             label = length(df[df$type == "Different TU" & df$int >= 200,]$int), size = 2.5)
  return(plt)
}
type_vs_interface(gene_pairs)

int_perc <- function(df, i, greater_than=FALSE){
  # Calc perc of interface-forming pairs for a given intervening gene count.
  if (greater_than == FALSE){
    nom <- length(filter(df, intervening == i, int >= 200)$int)
    denom <- length(filter(df, intervening == i, int < 200)$int) + nom
  } else {
    nom <- length(filter(df, intervening >= i, int >= 200)$int)
    denom <- length(filter(df, intervening >= i, int < 200)$int) + nom
  }
  return(list(nom, denom, nom/denom*100)) # last element: % forming interface
}

get_errors <- function(df){
  # Get binomial confidence intervals for data frame
  yd <- list(int_perc(df, 0), int_perc(df, 1), 
             int_perc(df, 2), int_perc(df, 3, TRUE))
  yvals <- as.numeric(c(yd[[1]][3], yd[[2]][3], yd[[3]][3], yd[[4]][3]))      
  xvals <- c('0', '1', '2', '>2')
  counts <- as.numeric(c(yd[[1]][2], yd[[2]][2], yd[[3]][2], yd[[4]][2]))
  yhi<- c()
  ylo <- c()
  for (i in yd){
    ylo <- append(ylo, 
                  binom.confint(as.numeric(i[1]), 
                                as.numeric(i[2]), 0.68)[[5]][11])
    yhi <- append(yhi, 
                  binom.confint(as.numeric(i[1]), 
                                as.numeric(i[2]), 0.68)[[6]][11])
  }
  df <- data.frame(xvals, yvals, yhi, ylo, counts)
  df$xvals <- factor(df$xvals, levels = xvals)
  df$type <- c('Adjacent', 'Non-adjacent',
               'Non-adjacent', 'Non-adjacent')
  return(df)
}

barplot_intervening_interface <- function(df){
  # Save at 3.5x3.5 inches
  plt_data <- get_errors(df)
  plt <- ggplot(plt_data, aes(xvals, yvals, fill=type)) +
    geom_bar(stat='identity', colour = 'black', lwd = 0.4) + 
    geom_errorbar(aes(ymax = yhi*100, ymin=ylo*100),
                  width=0.1) +
    scale_fill_manual(values = c('darkorange', 'skyblue3')) +
    xlab('Intervening genes between\ninteracting pair') +
    ylab('Gene pairs forming\nan interface (%)') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    theme(text = element_text(size = 8),
          legend.key.size = unit(0.25, "cm"),
          legend.justification = 'right', 
          legend.position=c(1,0.85)) +
    annotate("text", x = 1, y = 7, label = plt_data$count[1], size = 2.5) +
    annotate("text", x = 2, y = 7, label = plt_data$count[2], size = 2.5) +
    annotate("text", x = 3, y = 7, label = plt_data$count[3], size = 2.5) +
    annotate("text", x = 4, y = 7, label = plt_data$count[4], size = 2.5) +
    guides(fill=guide_legend(title="Gene pairs"))
  return(plt)
}

barplot_intervening_interface(comp1)
barplot_intervening_interface(ops_only)
type_vs_interface(gene_pairs)
#############################################################################


boxplot_intervening_interface <- function(df){
  # Returns boxplot of interface size vs pos. of gene in or out of operon.
  df <- mutate(df, int_facs = TRUE)
  df[df$type == 'Adjacent' & df$intervening == 0,]$int_facs <- '0'
  df[df$type == 'Non-adjacent' & df$intervening == 1,]$int_facs <- '1'
  df[df$type == 'Non-adjacent' & df$intervening > 1,]$int_facs <- '>1'
  df[df$type ==  "Different TU",]$int_facs <- "Different TU"
  df$int_facs <- factor(df$int_facs, levels = c('0', '1', '>1', "Different TU"))
  plt <- ggplot(df[df$int >= 200,], aes(int_facs, int)) +  
    geom_boxplot(fill = c('darkorange', 'skyblue3', 'skyblue3', '#C22326'), 
                 outlier.size = 1) +
    xlab('Intervening genes\n') + ylab(expression(paste("Interface size (\uc5"^"2",")"))) +
    theme(text = element_text(size=8), 
          axis.title.x = element_text(vjust = -0.5)) +
    scale_y_log10(breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000)) +
    annotate("text", x = 1.2, y = 500, 
             label = length(df[df$int_facs == "0" & df$int >= 200,]$int), size = 2.5) +
    annotate("text", x = 2.2, y = 500, 
             label = length(df[df$int_facs == "1" & df$int >= 200,]$int), size = 2.5) +
    annotate("text", x = 3.2, y = 500, 
             label = length(df[df$int_facs == ">1" & df$int >= 200,]$int), size = 2.5) +
    annotate("text", x = 4.2, y = 500, 
             label = length(df[df$int_facs == "Different TU" & df$int >= 200,]$int), size = 2.5)
  return(plt)
}

bxoplot_intervening_interface(gene_pairs)


#############################################################################
# E.coli Y2H PPI figures

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

#############################################################################
# Interface vs intervening, stacked barplot

stackedPlotter <- function(df){
  xlevs <- c("0", "1", "2", "3", "4", "5", "6", "", ">6")
  df <- mutate(df, int_existence = "Yes", col = "darkorange")
  df[df$int < 200,]$int_existence <- "No"
  df[df$int_existence == "No",]$col <- "darkorange2"
  df$int_existence <- factor(df$int_existence, levels = c("Yes", "No"))
  df[df$intervening > 6,]$intervening <- ">6"
  df$intervening <- factor(df$intervening, levels = xlevs)
  plt <- ggplot() +
    geom_bar(data = filter(df, intervening == 0), aes(intervening, fill = col)) +
    geom_bar(data = filter(df, intervening != 0), aes(intervening, fill = int_existence)) +
    scale_fill_manual(values = c("darkorange1", "goldenrod1", "skyblue1", "deepskyblue3")) +
    scale_x_discrete(limits = xlevs) +
    theme(legend.position = "none",
          text = element_text(size = 8)) +
    xlab('Intervening genes between pair') +
    ylab('Gene pairs') +
    geom_vline(xintercept = c(8), linetype = 2, alpha = 0.5, lwd = 0.6)
  return(plt)
}
stackedPlotter(ops_only)
