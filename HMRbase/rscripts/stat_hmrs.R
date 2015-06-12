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



dir_all <- "~/data/hmr/all"
file_out <- "~/data/hmr/all/statistics.txt"
files <- list.files(path = dir_all, pattern = "*.hmr", all.files = FALSE,
                    full.names = FALSE, recursive = TRUE,
                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

num_files <- length(files)

stat <- data.frame(name=files, row.names = files)
stat$length_total <- 0
stat$length_mean <- 0
for (fnm in files) {
  f <- file.path(dir_all, fnm)
  hmr <- read.table(f, header=FALSE, sep='\t')
  hmr$size <- (hmr[, 3] - hmr[, 2]+ matrix(1, nrow=dim(hmr)[1], ncol=1))
  stat[fnm, 'length_total'] <- sum(hmr$size)
  stat[fnm, 'length_mean'] <- mean(hmr$size)
}
stat$tissue <- sapply(rownames(stat), function(x) strsplit(x, '/')[[1]][1])
stat$sample <- sapply(rownames(stat), function(x) strsplit(x, '/')[[1]][2])
stat$sample <- sapply(stat$sample, function(x) strsplit(x, '.hmr')[[1]][1])
write.table(stat, file=file_out, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

tissuecol <- data.frame(tissue=unique(stat$tissue), row.names=unique(stat$tissue))
tissuecol$color <- rainbow(dim(tissuecol)[1])

p1 <- ggplot(stat, aes(sample, length_total))
p1 <- p1 + geom_bar(width=.3, fill=tissuecol[stat$tissue, 'color']) + coord_flip()
p1 <- p1 + theme(axis.text.y = element_text(colour='black', size=9))
p2 <- ggplot(stat, aes(sample, length_mean))
p2 <- p2 + geom_bar(width=.3, fill=tissuecol[stat$tissue, 'color']) + coord_flip()
p2 <- p2 + theme(axis.text.y = element_text(colour='black', size=9))
pdf(file=paste0(file_out, '-fig.pdf'), width=10, height=10, pointsize=10)
multiplot(p1, p2, cols=2)
dev.off()