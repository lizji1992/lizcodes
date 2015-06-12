library(plyr)
library(mclust)

###############################################################################
dir.signal <- '~/codes/lizcodes/subHMR/test/meth'
dir.out <- '~/data/subHMR/cpg_cluster'
file.index <- '~/codes/lizcodes/subHMR/test/cpgs.bed'
samplenum <- .Machine$integer.max
###############################################################################

source('~/codes/HMRbase/rscripts/lib/MethFiles.R')
cpgtab <- loadIndex(file.index)

infiles <- list.files(path=dir.signal)
for (f in infiles) {
  sampleid <- strsplit(f, '\\.')[[1]][1]
  signaltab <- loadMeth(file.path(dir.signal, f), sampleid)
  cpgtab <- merge(cpgtab, signaltab, by.x=c('chr', 'start'),
                all.x=TRUE, all.y=FALSE)
}

cpgtab <- cpgtab[sapply(1:dim(cpgtab)[1], function(x) all(!is.na(cpgtab[x, 3:(2+length(infiles))]))), ]
d <- dist(cpgtab, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit)