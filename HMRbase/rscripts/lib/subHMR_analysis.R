

assign.group <- function(tab) {
  tab.op <- tab[, c('chr', 'start', 'end')]
  
  num <- dim(tab.op)[1]
  
  tab.op$lastend <- 0
  tab.op[2:num, 'lastend'] <- tab.op[1:num-1, 'end']
  tab.op$gap <- (tab.op$start - tab.op$lastend)
  
  collapsed_groups <- rle(tab.op$gap)$lengths
  tab$group <- unlist(sapply(1:length(collapsed_groups), function(x) rep(x, collapsed_groups[x])))
  tab$num_member <- unlist(sapply(1:length(collapsed_groups), function(x) rep(collapsed_groups[x], collapsed_groups[x])))
  
  return(tab)
}