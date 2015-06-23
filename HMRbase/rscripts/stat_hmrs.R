library(plyr)
library(ggplot2)
library(grid)
library(RColorBrewer)

############
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
############

get_stat <- function(dir_all, fire_out) {
  files <- list.files(path = dir_all, pattern = "*.hmr", all.files = FALSE,
                      full.names = FALSE, recursive = TRUE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  num_files <- length(files)
  
  stat <- data.frame(name=files, row.names = files)
  stat$length_total <- 0
  stat$length_mean <- 0
  stat$num <- 0
  for (fnm in files) {
    f <- file.path(dir_all, fnm)
    hmr <- read.table(f, header=FALSE, sep='\t')
    hmr$size <- (hmr[, 3] - hmr[, 2]+ matrix(1, nrow=dim(hmr)[1], ncol=1))
    stat[fnm, 'num'] <- dim(hmr)[1]
    stat[fnm, 'length_total'] <- sum(hmr$size)
    stat[fnm, 'length_mean'] <- mean(hmr$size)
  }
  stat$tissue <- sapply(rownames(stat), function(x) strsplit(x, '_')[[1]][1])
  stat$sample <- sapply(rownames(stat), function(x) strsplit(x, '.hmr')[[1]][1])
  write.table(stat, file=file_out, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
  return(stat)
}


list.sizes <- function(dir_all) {
  files <- list.files(path = dir_all, pattern = "*.hmr", all.files = FALSE,
                      full.names = FALSE, recursive = TRUE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  num_files <- length(files)
  
  l <- matrix(, nrow = 0, ncol = 1)
  for (fnm in files) {
    f <- file.path(dir_all, fnm)
    hmr <- read.table(f, header=FALSE, sep='\t')
    hmr$size <- (hmr[, 3] - hmr[, 2]+ matrix(1, nrow=dim(hmr)[1], ncol=1))
    l <- rbind(l, hmr$size) 
  }
  return(l)
}

############

stat1 <- get_stat("~/data/hmr/choose", "~/data/hmr/choose/statistics.txt")
stat2 <- get_stat("~/data/hmr/choose_ct", "~/data/hmr/choose_ct/statistics.txt")

sizes1 <- list.sizes("~/data/hmr/choose")
sizes2 <- list.sizes("~/data/hmr/choose_ct")
stat <- merge(stat1[, c('num', 'length_mean', 'tissue', 'sample')], stat2[, c('num', 'length_mean')],
              by="row.names")

sizes <- data.frame(size=rbind(sizes1, sizes2), row.names = NULL)
sizes$model <- c(rep('hmm', dim(sizes1)[1]), rep('cthmm', dim(sizes2)[1]))


stat.num <- rbind(stat1[, c('num', 'tissue', 'sample')], stat2[, c('num', 'tissue', 'sample')])
stat.num$model <- c(rep('hmm', dim(stat1)[1]), rep('cthmm', dim(stat1)[1]))


p1 <- ggplot(data=stat.num, aes(stat.num$sample, stat.num$num, fill=stat.num$model))
p1 <- p1 + geom_bar(width=.3, stat="identity", position="dodge") + coord_flip() + theme_bw()
p1 <- p1 + theme(axis.title.y=element_blank(), legend.title=element_blank())
p1

stat.length_mean <- rbind(stat1[, c('length_mean', 'tissue', 'sample')], stat2[, c('length_mean', 'tissue', 'sample')])
stat.length_mean$model <- c(rep('hmm', dim(stat1)[1]), rep('cthmm', dim(stat1)[1]))
p2 <- ggplot(data=stat.length_mean, aes(stat.length_mean$sample, stat.length_mean$length_mean, fill=stat.length_mean$model))
p2 <- p2 + geom_bar(width=.3, stat="identity", position="dodge") + coord_flip() + theme_bw()
p2 <- p2 + theme(axis.title.y=element_blank(), legend.title=element_blank())
p2

pdf(file=paste0(file_out, '-fig.pdf'), width=10, height=10, pointsize=10)
multiplot(p1, p2, cols=2)
dev.off()

p <- ggplot(data=sizes, aes(sizes$size, y=..count.., fill=sizes$model)) +
  geom_density(alpha=0.4) + theme_bw()
p

pdf(file=paste0(file_out, '-size_distr20000-70000.pdf'), width=10, height=10, pointsize=10)
p + scale_x_continuous(limits=c(20000, 70000))
dev.off()

