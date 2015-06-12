require(graphics)
require(grDevices)
library(gplots)
library(RColorBrewer)

dropbox <- '~/data'
f <- file.path(dropbox, 'subhmr/f0.125-s10-d1/f0.125-s10-d1-freq')

tab <- read.table(paste0(f, '.tab'), skip=0, sep=' ', header = TRUE)

tab.sorted <- tab[with(tab, order(-freq)),]
write.table(tab.sorted,
            file =paste0(f, '.id.sorted'),
            quote=FALSE, sep='\t', row.names=FALSE)
###########################


range <- 0.8
top <- 500
###########################


tab.range <- tab.sorted[(tab.sorted[, 'freq'] < range), ]
tab.top <- tab.range[1:top, 3:dim(tab.range)[2]]
rownames(tab.top) <- tab.top$id


############################################

mat  <- as.matrix(tab.top)
col.group <- brewer.pal(5,"Blues")

pdf(file=paste0(f, "_r", range, "_t", top, ".heatmap.pdf"),width=15, height=10, pointsize=10)
layout(t(1:1))
heatmap.2(mat, dendrogram = c("column"), scale = "row", trace = "none", na.rm = TRUE,
          col = c("cornsilk", "chocolate1"),
          labRow = NA, margin=c(15, 15),
          #ColSideColors = c( rep(col.group[1], 1), rep(col.group[2], 4), rep(col.group[3], 3),
          #                  rep(col.group[4], 9), rep(col.group[5], 16) ),
          keysize = 0.5, key.title = "occurence in samples", key.xlab = "1 or 0", key.ylab=NA,
          main = paste0("cluster samples by cell type"),
          xlab = "Samples", ylab = "subHMRs")
par(lend = 1)
legend("topright", c("Sperm", "iPSC", "ESC", "Blood", "Brain"), 
       col = col.group, lty = 1, lwd = 10)
dev.off()