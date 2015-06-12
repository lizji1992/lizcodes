library(plyr)

loadIndex <- function(fin) {
  tab <- read.table(fin, header=FALSE, sep='\t',
                    colClasses=c('character', 'numeric'), col.names=c('chr', 'start'))
  return(tab)
}

loadMeth <- function(fin, sampleid) {
  tab <- read.table(fin, header=FALSE, sep='\t', row.names=NULL,
                    colClasses=c('character', 'numeric', rep('NULL', 2), rep('numeric', 2)),
                    col.names=c('chr', 'start', rep('', 2), sampleid, 'reads_coverage'))
  tab <- tab[tab[, 'reads_coverage'] > 0, c('chr', 'start', sampleid)]
  return(tab)
}
