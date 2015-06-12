library(mclust)

dropbox <- '~/data'
f <- file.path(dropbox, 'subHMR/panels/panels-100.signal')
outpref <- file.path(dropbox, 'subHMR/panels/Figure/panels-1000samples/panels-1000samples')
signal_all <- read.table(paste0(f), header=FALSE, sep=',')

num_sample <- 1000
num_group <- 5
signal <- signal_all[sample(nrow(signal_all), num_sample), ]
x <- 1:100

d <- dist(signal, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
pdf(file=paste0(outpref, '_tree.pdf'), width=10, height=6, pointsize=10)
layout(t(1:1))
plot(fit);
dev.off()

groups <- cutree(fit, h=5)
num <- length(unique(groups[]))
for (k in 1:num) {
  fnm <- paste(outpref, num_group, k, sep = '-')
  idx <- as.numeric(names(groups[groups[] == k]))
  group <- signal_all[idx,]
  pdf(file=paste0(fnm, '.pdf'), width=10, height=6, pointsize=10)
  layout(t(1:1))
  group.mean <- colMeans(group)
  plot(x, group.mean, col='red', xlab="position", ylab="signal", type="l")
  title(main=paste0(length(idx), " islands among ", num_sample, " in total." ))
  dev.off()
}



signal.mean <- colMeans(signal)
plot(x, signal.mean)
y <- as.matrix(signal[500,])
plot( 1:100, y[1,], col='red', xlab="position", ylab="signal", type="l")

plot( x, signal[25,], col='red', xlab="position", ylab="signal", type="p")

# Ward Hierarchical Clustering
fit <- hclust(d, method="ward")
group1 <- signal[as.numeric(names(groups[groups[] == 88])),]
group1.mean <- colMeans(group1)
plot( x, group1.mean, col='red', xlab="position", ylab="signal", type="l")




plot( x, group1[11,], col='red', xlab="position", ylab="signal", type="l")


# mclust
fit.mclust <- Mclust(signal)


plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")


heatmap.2(as.matrix(signal), dendrogram = c("column"), scale = "row", trace = "none", na.rm = TRUE,
          col = c("cornsilk", "chocolate1"),
          labRow = NA, margin=c(15, 15),      
          keysize = 0.5, key.title = "occurence in samples", key.xlab = "1 or 0", key.ylab=NA,
          main = paste0("cluster samples by cell type"),
          xlab = "Samples", ylab = "subHMRs")

plot( 1:100, signal.mean, col='blue', xlab="position", ylab="signal", type="p")
lo <- loess(y[1,] ~ x)

plot(x, y)
lines(predict(lo), col='red', lwd=2)