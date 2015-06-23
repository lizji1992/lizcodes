library(plyr)

loadSHMR <- function(tabfile) {
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c('character',rep('numeric',2), 'character', 'numeric', 'NULL'),
                    col.names=c('chr', 'start', 'end', 'state', 'freq', ''))
  tab$size <- (tab$end - tab$start + matrix(1, nrow=dim(tab)[1], ncol=1))
  return(tab)
}

loadBED_chromHMM <- function(tabfile) {
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c('character',rep('numeric',2), 'character', 'numeric'),
                    col.names=c('chr', 'start', 'end', 'state', 'num_cpg'))
  tab$size <- (tab$end - tab$start + matrix(1, nrow=dim(tab)[1], ncol=1))
  return(tab)
}


loadSHMR <- function(tabfile, list.sample, osep=',') {
  tabfile <- file.tab
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c('character', rep('numeric', 4), 'character', rep('numeric', 2)),
                    col.names=c('chr', 'start', 'end', 'istart', 'iend', 'state0', 'num_cpg', 'score'))
  tab$size <- (tab$end - tab$start + matrix(1, nrow=dim(tab)[1], ncol=1))
  
  oldstates <- unique(tab$state0)
  split_state <- strsplit(oldstates, osep)
  state <- data.frame(matrix( as.numeric(unlist(split_state)), nrow=length(split_state), byrow=T),
                      row.names=oldstates, stringsAsFactors=FALSE)
  colnames(state) <- 1:length(state)
  
  idx.shuffle <- list.sample[, c('index0', 'index')]
  rownames(idx.shuffle) <- idx.shuffle$index
  state <- data.frame(sapply(1:length(state), function(x) state[, idx.shuffle[x, 'index0']]))
  state$newstate <- sapply(1:dim(state)[1], function(x) paste0(state[x, ], collapse=''))
  rownames(state) <- oldstates
  tab$state <- state[tab$state0, 'newstate']
  return(tab)
}


loadSHMR_tss <- function(tssfile, pfile, list.sample, osep=',') {
  tab <- read.table(tssfile, header=FALSE, sep='\t',
                    colClasses=c('character', rep('numeric', 4), 'character', 'numeric', 'NULL',
                                 'NULL', 'numeric', rep('NULL', 3), 'character'),
                    col.names=c('chr', 'start', 'end', 'istart', 'iend', 'state0', 'num_cpg', '',
                                '', 'tss', rep('', 3), 'strand'))  
  tab <- tab[!duplicated(tab[, c('start')]), ]
  rownames(tab) <- NULL
  tab.p <- read.table(pfile, header=FALSE, sep='\t',
                    colClasses=c(rep('NULL', 6), 
                                 'character', rep('numeric', 4), rep('NULL', 3)),
                    col.names=c(rep('', 6), 
                                'chr', 'start', 'end', 'istart', 'iend', rep('', 3)))
  tab.p <- tab.p[!duplicated(tab.p[, c('start')]), ]
  
  tab.p$promoter <- TRUE
  tab <- join(tab, tab.p, by=c('chr', 'start', 'end', 'istart', 'iend'))
  tab[is.na(tab[, 'promoter']), 'promoter'] <- FALSE

  
  tab$size <- (tab$end - tab$start + matrix(1, nrow=dim(tab)[1], ncol=1))
  oldstates <- unique(tab$state0)
  split_state <- strsplit(oldstates, osep)
  state <- data.frame(matrix( as.numeric(unlist(split_state)), nrow=length(split_state), byrow=T),
                      row.names=oldstates, stringsAsFactors=FALSE)
  colnames(state) <- 1:length(state)
  
  idx.shuffle <- list.sample[, c('index0', 'index')]
  rownames(idx.shuffle) <- idx.shuffle$index
  state <- data.frame(sapply(1:length(state), function(x) state[, idx.shuffle[x, 'index0']]))
  state$newstate <- sapply(1:dim(state)[1], function(x) paste0(state[x, ], collapse=''))
  rownames(state) <- oldstates
  tab$state <- state[tab$state0, 'newstate']
  rownames(tab) <- NULL
  return(tab)
}


loadSHMR_p <- function(tabfile, tssfile, list.sample, osep=',') {
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c(rep('NULL', 6), 
                                 'character', rep('numeric', 4), 'character', 'numeric', 'NULL'),
                    col.names=c(rep('', 6), 
                                'chr', 'start', 'end', 'istart', 'iend', 'state0', 'num_cpg', ''))
  tab <- tab[!duplicated(tab[, c('start')]), ]
  
  oldstates <- unique(tab$state0)
  split_state <- strsplit(oldstates, osep)
  state <- data.frame(matrix( as.numeric(unlist(split_state)), nrow=length(split_state), byrow=T),
                      row.names=oldstates, stringsAsFactors=FALSE)
  colnames(state) <- 1:length(state)
  
  idx.shuffle <- list.sample[, c('index0', 'index')]
  rownames(idx.shuffle) <- idx.shuffle$index
  state <- data.frame(sapply(1:length(state), function(x) state[, idx.shuffle[x, 'index0']]))
  state$newstate <- sapply(1:dim(state)[1], function(x) paste0(state[x, ], collapse=''))
  rownames(state) <- oldstates
  tab$state <- state[tab$state0, 'newstate']
  rownames(tab) <- NULL
  return(tab)
}

loadSHMR_np <- function(tabfile, list.sample, osep=',') {
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c(rep('NULL', 6), 
                                 'character', rep('numeric', 4), 'character', 'numeric', 'NULL'),
                    col.names=c(rep('', 6), 
                                'chr', 'start', 'end', 'pstart', 'pend', 'state0', 'num_cpg', ''))
  tab <- tab[!duplicated(tab[, c('start')]), ]
  tab$size <- (tab$end - tab$start + matrix(1, nrow=dim(tab)[1], ncol=1))
  
  oldstates <- unique(tab$state0)
  split_state <- strsplit(oldstates, osep)
  state <- data.frame(matrix( as.numeric(unlist(split_state)), nrow=length(split_state), byrow=T),
                      row.names=oldstates, stringsAsFactors=FALSE)
  colnames(state) <- 1:length(state)
  
  idx.shuffle <- list.sample[, c('index0', 'index')]
  rownames(idx.shuffle) <- idx.shuffle$index
  state <- data.frame(sapply(1:length(state), function(x) state[, idx.shuffle[x, 'index0']]))
  state$newstate <- sapply(1:dim(state)[1], function(x) paste0(state[x, ], collapse=''))
  rownames(state) <- oldstates
  tab$state <- state[tab$state0, 'newstate']
  rownames(tab.p) <- NULL
  return(tab)
}


loadSampleList <- function(samplelist_file, sample_id_file, category_file='') {

  list.all <- read.table(samplelist_file, header=FALSE, sep='\t', col.names=c('proj', 'sample_id'))
  list.all$proj <- sapply(list.all$proj, function(x) tail(strsplit(as.character(x), "/")[[1]], n=1))
  list.all$proj <- sapply(list.all$proj, function(x) strsplit(as.character(x), "\\.")[[1]][1])
  list.all$proj <- sapply(list.all$proj, function(x) strsplit(as.character(x), "_")[[1]][2])
  row.names(list.all) <- list.all$sample_id
  
  map.category <- data.frame(sample_id=as.character(list.all$sample_id), cat_id=as.character(list.all$sample_id),
                             stringsAsFactors=FALSE, row.names=list.all$sample_id) 
  if (category_file != '') {
    list.category <- read.table(category_file, header=FALSE, sep='\t', col.names=c('cat_id', 'sample_id'),
                                stringsAsFactors=FALSE, colClasses=c(rep('character', 2)))
    rownames(list.category) <- list.category$sample_id
    map.category[rownames(list.category), 'cat_id'] <- as.character(list.category$cat_id)
  }
  
  list.sample <- read.table(sample_id_file, header=FALSE, sep='\t', col.names=c('index0', 'sample_id'), stringsAsFactors=FALSE)
  list.sample$index0 <- (list.sample$index0 + matrix(1, nrow=dim(list.sample)[1], ncol=1))
  row.names(list.sample) <- list.sample$index0
  list.sample$project <- list.all[list.sample$sample_id, 'proj']
  list.sample$sample_id <- map.category[list.sample$sample_id, 'cat_id']
  list.sample$tissue <- sapply(list.sample$sample, function(x) head(strsplit(as.character(x), "_")[[1]], n=1))
  list.tissue <- data.frame(unique(list.sample$tissue))
  list.tissue$id <- rownames(list.tissue)
  row.names(list.tissue) <- list.tissue[,1]
  list.sample$tissue_id <- list.tissue[list.sample$tissue, 'id']
  list.sample <- list.sample[order(list.sample$tissue_id, list.sample$sample_id), ]
  list.sample$index <- 1:dim(list.sample)[1]
  rownames(list.sample) <- list.sample$sample_id 
  return(list.sample)
}

loadIntervals_UCSC <- function(infile) {
  list.gene <- read.table(infile, header=TRUE, sep='\t', colClasses=c(rep('character', 3) , 'NULL'))
  list.gene$start <- sapply(list.gene[, 3], function(x) head(strsplit(as.character(x), "-")[[1]], n=1))
  list.gene$start <- sapply(list.gene$start, function(x) tail(strsplit(as.character(x), ":")[[1]], n=1))
  list.gene$end <- sapply(list.gene[, 3], function(x) tail(strsplit(as.character(x), "-")[[1]], n=1))
  list.gene$start <- sapply(list.gene$start, function(x) as.numeric(gsub(",", "", x)))
  list.gene$end <- sapply(list.gene$end, function(x) as.numeric(gsub(",", "", x)))
  row.names(list.gene) <- list.gene$gene
  list.gene[, 3] <- NULL
  return(list.gene)
}  