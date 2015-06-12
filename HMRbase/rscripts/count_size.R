require(plotrix)
require(graphics)

x1 <- read.table('~/data/subhmr/f0.125-s10-d1/f0.125-s10-d1-freq', header=FALSE, sep='\t')
pdf(file="~/data/subhmr/f0.125-s10-d1/f0.125-s10-d1-freq_threshold.pdf",width=10, height=6, pointsize=10)


x1$size <- x1[, 3] - x1[, 2]


l <- list (x1$size)

layout(t(1:1))
multhist(l, freq=FALSE, breaks=400, xlim=c(0, 100),
         xlab="subHMR size (bp)", ylab="density", main="size of subHMRs under different thresholds")
legend("topright", legend=c(0.9, 0.85, 0.8, 0.7, 0.6), fill=gray.colors(5),
       text.font=0.5)
dev.off()
#hist(x$size, freq=FALSE, breaks=500, xlim=c(0, 10000), plot=TRU)


x6=read.table('~/Dropbox/subHMRs/hmr/Blood_1.hmr', header=FALSE, sep='\t')
x6$size <- x6[, 3] - x6[, 2]
hist(x6$size, freq=FALSE, breaks=400,
         xlab="HMR size (bp)", ylab="density", main="size of HMR")