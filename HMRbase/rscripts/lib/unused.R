
#for (i in 1:num_sample) {
#  tab.select <- tab[tab[,"count"]==i, ]
#  state <- unique(tab.select$state)
#  stat <- data.frame(count=matrix(numeric(), nrow=length(state)), row.names=state)
#  stat <- stat.count(tab.select, stat)
#  stat <- stat.size(tab.select, stat)
#  stat <- stat[order(-stat[, "len_mean"]),]
#  barplot(stat[, 'len_mean'], main=c("Mean length (density:", i, ")"))
#  axis(1, at = 1:length(state), labels=rownames(stat), las=2, cex.axis=0.8)
#}


