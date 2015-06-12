library(plyr)
library(Matrix)

state.merge <- function(statelist, merge) {
  statelist.new <- vector(mode="character", length=length(merge))
  for (i in 1:length(merge)) {
    if ('1' %in% statelist[merge[[i]]]) {
      statelist.new[i] <- '1'
    }
  }
  return(statelist.new)
}

state.assign_tissue_id <- function(state, list.sample, osep='', tsep=',') {
  state.uniq <- data.frame(unique(state))
  num_state <- length(state.uniq)
  rownames(state.uniq) <- state.uniq[, 1]
  num_sample <- length(strsplit(as.character(state.uniq[1,1]), osep)[[1]])
  for (i in 1:num_sample) {
    state.uniq[, i] <- '0'
  }
  colnames(state.uniq) <- list.sample[1:num_sample, 'index']
  split_state <- strsplit(rownames(state.uniq), osep)
  state.uniq[, ] <- data.frame(matrix( as.numeric(unlist(split_state)), nrow=length(split_state), byrow=T ))  
  
  for (i in 1:num_sample) {
    state.uniq[state.uniq[, i] != 0, i] <- list.sample[i, 'tissue_id']
  }
  
  offset <- 0
  tissue <- unique(list.sample$tissue_id)
  for (t in tissue) {
    if (offset < length(tissue)-1) {
      idx <- list.sample[list.sample[, 'tissue_id']==t, 'index']
      state.uniq <- cbind(state.uniq[, 1:(idx[length(idx)]+offset)],
                          data.frame(matrix(rep(tsep, num_state), nrow=num_state), stringsAsFactors=FALSE),
                          state.uniq[, (idx[length(idx)]+1+offset):length(state.uniq)] )
      offset <- offset + 1;
    }
  }
  state.uniq$state <- sapply(rownames(state.uniq), function(x) paste0(state.uniq[x, ], collapse=''))
  state.uniq$oldstate <- rownames(state.uniq)
  return(state.uniq[state, c('oldstate', 'state')])
}

#==============================================================================
pos.stateorder <- function(tab, stat, item, midstate, range=20) {
  states <- rownames(stat)
  posmat <- data.frame(matrix(0, nrow=length(states), ncol=range*2), row.names = rownames(stat))
  colnames(posmat) <- c(-range:-1, 1:range)
  
  num <- dim(tab)[1]
  tab$idx <- 1:num
  mididx <- tab[tab[, 'state']==midstate, 'idx']
  for (i in 1:length(mididx)) {
    notover <- TRUE
    notover_left <- TRUE
    notover_right <-TRUE
    offset <- 1
    mid <- mididx[i]
    strd <- 1
    if (tab[mid, 'strand'] == '-') {
      strd <- -1
    }
    while(notover) {
      if (offset <= range && (mid - offset) > 0
          && tab[(mid - offset), 'istart'] == tab[mid, 'istart']
          && (mid-offset) > mididx[max(1, i-1)] ) {
        state <- tab[(mid - offset), 'state']
        posmat[state, as.character(-strd*offset)] <- posmat[state, as.character(-strd*offset)] + tab[(mid - offset), item]
      } else {
        notover_left <- FALSE
      }
      if (offset <= range && (mid + offset) <= num
          && tab[(mid + offset), 'istart'] == tab[mid, 'istart']
          && (mid+offset) < mididx[min(i+1, length(mididx))] ) {
        state <- tab[(mid + offset), 'state']
        posmat[state, as.character(strd*offset)] <- posmat[state, as.character(strd*offset)] + tab[(mid + offset), item]
      } else {
        notover_right <- FALSE
      }
      notover <- notover_left | notover_right
      offset <- offset + 1
    }
  }
  posmat <- posmat[, colSums(posmat) > 0]
  posmat <- posmat[setdiff(rownames(posmat), midstate), ]
  order <- data.frame(matrix('', nrow=1, ncol=dim(posmat)[2]), row.names = c("state"), stringsAsFactors=FALSE)
  colnames(order) <- colnames(posmat)
  order['state', ] <- sapply(1:dim(posmat)[2], function(x) rownames(posmat)[which.max( posmat[, x] )])
  return(list(posmat, order))
}


strip_pos.sortitem <- function(stat, item) {
  state <- rownames(stat)
  strip_pos <- data.frame(val=stat[, item], state=state, row.names = state)
  strip_pos <- strip_pos[order(strip_pos$val), ]
  strip_pos$idx <- 1:dim(strip_pos)[1]
  return(strip_pos[rownames(stat), 'idx'])
}

strip_pos.state <- function(state, osep='') {
  state.uniq <- unique(state)
  strip_pos <- data.frame(state.uniq)
  rownames(strip_pos) <- strip_pos[,1]
  num_state <- length(state.uniq)
  num_sample <- length(strsplit(as.character(state[1]), osep)[[1]])
  for (i in 1:num_sample) {
    strip_pos[, i] <- '0'
  }
  split_state <- strsplit(rownames(strip_pos), osep)
  strip_pos[, ] <- data.frame(matrix( as.numeric(unlist(split_state)), nrow=length(split_state), byrow=T ))
  colnames(strip_pos) <- 1:num_sample
  strip_pos$freq <- sapply(1:num_state, function(x) nnzero(strip_pos[x, 1:num_sample ]) )
  
  mid <- num_sample %/% 2
  list.lscore <- data.frame(matrix(0, nrow=num_sample, ncol=1), row.names=1:num_sample)
  list.lscore[1:mid, 1] <- 1:mid
  list.rscore <- data.frame(matrix(0, nrow=num_sample, ncol=1), row.names=1:num_sample)
  list.rscore[(num_sample-mid+1):num_sample, 1] <- mid:1

  strip_pos$lscore <- sapply(1:num_state, function(x) sum(list.lscore[ strip_pos[x, 1:num_sample ]!=0, 1] ) )
  strip_pos$rscore <- sapply(1:num_state, function(x) sum(list.rscore[ strip_pos[x, 1:num_sample]!=0, 1] ) )
  strip_pos$side <- 'l'
  for (i in 1:num_state) {
    if (strip_pos[i, 'rscore'] > strip_pos[i, 'lscore']) {
      strip_pos[i, 'side'] <- 'r'
    }
    else if (strip_pos[i, 'rscore'] == strip_pos[i, 'lscore']) {
      strip_pos[i, 'side'] <- sample(c('l','r'), 1)
    }
  }
  fullstate <- paste(rep('1', num_sample), collapse=',')
  if (fullstate %in% rownames(state.uniq)) {
    strip_pos[fullstate, 'side'] <- 'm'
  }
  strip_pos.side <- lapply(split(strip_pos, strip_pos$side), function(x) data.frame(x))
  #strip_pos.side[[1]]$score <- (strip_pos.side[[1]]$lscore * strip_pos.side[[1]]$freq - strip_pos.side[[1]]$rscore)
  for (i in 1:length(strip_pos.side)) {
    if (strip_pos.side[[i]][1, 'side'] == 'l') {
      strip_pos.side[[i]] <- strip_pos.side[[i]][order(strip_pos.side[[i]][, 'lscore'], strip_pos.side[[i]][, 'rscore']), ]
      strip_pos <- strip_pos.side[[i]]
    }
    else if (strip_pos.side[[i]][1, 'side'] == 'r') {
      strip_pos.side[[i]] <- strip_pos.side[[i]][order(-strip_pos.side[[i]][, 'rscore'], -strip_pos.side[[i]][, 'lscore']), ]
      strip_pos <- rbind(strip_pos, strip_pos.side[[i]])
    }
    else {
      strip_pos <- rbind(strip_pos, strip_pos.side[[i]])
    }
  }                                                                                   
  strip_pos$idx <- 1:num_state
  return(strip_pos[state.uniq, 'idx'])
}

#==============================================================================

hmrs.in_range <- function(tab, gene) {
  tab.gene <- tab[tab[,'chr']==gene$chr, ]
  tab.gene <- tab.gene[tab.gene[,'start'] >= gene$start, ]
  tab.gene <- tab.gene[tab.gene[,'end'] <= gene$end, ]
  return(tab.gene)
}

#==============================================================================
stat.count <- function(tab, stat) {
  num_fragment <- dim(tab)[1]
  state_count <- plyr::count(as.vector(tab), "state")
  row.names(state_count) <- state_count$state
  stat <- merge(stat, state_count, by="row.names") 
  rownames(stat) <- stat$Row.names
  stat$Row.names <- NULL
  return(stat)
}

stat.sum <- function(tab, stat, item) {
  for(s in rownames(stat)) {
    stat[s, paste0(item, '_total')] <- sum(tab[tab[, 'state']==s, item])
  }
  return(stat)
}

stat.mean <- function(tab, stat, item) {
  for(s in rownames(stat)) {
    stat[s, paste0(item, '_mean')] <- mean(tab[tab[, 'state']==s, item])
  }
  return(stat)
}

stat.size <- function(tab, stat) {
  for(s in rownames(stat)) {
    stat[s, 'len_total'] <- sum(tab[tab[,'state']==s, 'size'])
    stat[s, 'len_mean'] <- mean(tab[tab[,'state']==s, 'size'])
   # stat[s, 'len_sd'] <- sd(tab[tab[,'state']==s, 'size'])[1]
  }
  return(stat)
}
#==============================================================================



