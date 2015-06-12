library(RColorBrewer)
library(scales)
require(plotrix)
require(graphics)

dropbox <- '~/data'
f <- file.path(dropbox, 'subhmr/f0.125-s10-d1/f0.125-s10-d1-freq')

id <- read.table(paste0(f, '.id'), header=FALSE, sep='\t')
p <- read.table(paste0(f, '.p'), header=FALSE, sep='\t')
np <- read.table(paste0(f, '.np'), header=FALSE, sep='\t')
tss <- read.table(paste0(f, '.tssmid'), header=FALSE, sep='\t')


#################################################
p <- p[!duplicated(p[, 10]), ]
p.sizes <- p[, 3] - p[, 2]
p.sizes[!is.na(p.sizes)]

np <- np[!duplicated(np[, 4]), ]
np.sizes <- np[, 3] - np[, 2]
np.sizes[!is.na(np.sizes)]

pdf(file=paste0(f, ".sizes_pvsnp.pdf"),width=10, height=6, pointsize=10)
layout(t(1:1))


p.col <- rgb(1,0,0,0.2)
np.col <- rgb(0,0,1,0.2)

#xlim=range(min(c(p.sizes, np.sizes)), max(c(p.sizes, np.sizes)))

hist(p.sizes, col = p.col, xlim=c(0, 5000), main = 'Distribution of sizes of p-subHMRs and np-subHMRs', xlab = 'size (bp)')
hist(np.sizes, add=T, col = np.col, breaks=500)
## add a legend in the corner
legend('topright',c('p-subHMRs','np-subHMRs'), fill = c(p.col, np.col), bty = 'n', border = NA)

dev.off()

#################################################
tss <- tss[!duplicated(tss[, 4]), ]
p.id <- p[, 10]
np.id <- np[, 4]
p.tss <- tss[ tss[, 4] %in% p.id, ]
np.tss <- tss[ tss[, 4] %in% np.id, ]

pdf(file=paste0(f, ".TSS_pvsnp.pdf"),width=10, height=6, pointsize=10)
layout(matrix(c(1,2), 2, 1))

hist(p.tss[, 13], col='skyblue', xlim=c(0, 50000), main = 'Distribution of disatance to TSS of p-subHMRs',
     xlab = 'size (bp)', breaks=100)
hist(np.tss[, 13], xlim=c(0, 50000), col=scales::alpha('red',.5), main = 'Distribution of disatance to TSS of nonp-subHMRs',
     xlab = 'size (bp)', breaks=1000)
## add a legend in the corner
dev.off()

#################################################
pdf(file=paste0(f, ".freq_pvsnp.pdf"),width=10, height=6, pointsize=10)
layout(t(1:1))

densp <- density(p[, 11])
densnp <- density(np[, 5])

plot(densp$x, densp$y/max(densp$y), xlim=c(0,1), ylim=c(0,1.2), type='l', xlab = 'Enrichment',
     ylab = 'density(KDE)',
     main = 'Distribution of frequency of p-subHMRs and np-subHMRs', col='skyblue')
lines(densnp$x, densnp$y/max(densnp$y), xlim=c(0,1), col='red')
## add a legend in the corner
legend('topright',c('p-subHMRs','np-subHMRs'), fill = c('skyblue', 'red'), bty = 'n', border = NA)

dev.off()