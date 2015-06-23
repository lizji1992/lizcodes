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

loadBED_p <- function(tabfile) {
  p <- read.table(tabfile, header=FALSE, sep='\t',
                  colClasses=c('character',rep('numeric',2), 'character', rep('NULL', 5), 'character', 'numeric'),
                  col.names=c('chr', 'start', 'end', 'gene',  rep('', 5), 'id', 'num_cpg'))
  p <- p[!duplicated(p$id), ]
  return(p)
}

loadBED_np <- function(tabfile) {
  np <- read.table(tabfile, header=FALSE, sep='\t',
                   colClasses=c('character',rep('numeric',2), 'character', 'numeric'),
                   col.names=c('chr', 'start', 'end', 'id', 'num_cpg'))
  np <- np[!duplicated(np$id), ]  
  return(np)
}

loadBED_tss <- function(tabfile) {
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c('character', rep('numeric', 2), 'character', 'numeric', 'NULL',
                                 'NULL', 'numeric', rep('NULL', 2), 'character'),
                    col.names=c('chr', 'start', 'end', 'id', 'num_cpg', '',
                                '', 'tss', rep('', 2), 'strand'))
  tab <- tab[!duplicated(tab[, c('id')]), ]
  return(np)
}