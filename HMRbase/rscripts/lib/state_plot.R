library(gridExtra)
library(graphics)
library(plotrix)
library(ggplot2)
library(RColorBrewer)

library(plyr)
library(Matrix)

plothist.size.p <- function(stat.p, state.np) {
 # vals <- c(stat.p$size_total, stat.np$size_total)
#  vals <- vals[order(-vals)]
#  h1 <- vals[1]
#  h2 <- vals[1]*0.2
#  h3 <- vals[2]*1.01
  #stat.p <- stat.p[order(-stat.p[, size_total]), ]
  stat.p$sortstate <- factor(stat.p$state, levels =  stat.p$state[order(-stat.p$size_total)])
  stat.p$Context <- "Promoter"
  stat.p$color <- "Skyblue"
  stat.np$sortstate <- factor(stat.np$state, levels =  stat.np$state[order(-stat.np$size_total)])
  stat.np$Context <- "Non-promoter"
  stat.np$color <- "Red"
  
  color <- data.frame(color=c("Skyblue", "Red"), row.names=c("Promoter", "Non-promoter"))
  stat <- rbind(stat.p, stat.np)
  
  p <- ggplot(stat, aes(x = sortstate, y = size_total, fill=Context)) + 
    geom_bar(position="dodge") + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, colour='black'))
  p
  return(p)
}

plothist.cpg.p <- function(stat.p, state.np) {
  #stat.p <- stat.p[order(-stat.p[, num_cpg_total]), ]
  stat.p$sortstate <- factor(stat.p$state, levels =  stat.p$state[order(-stat.p$num_cpg_total)])
  stat.p$Context <- "Promoter"
  stat.p$color <- "Skyblue"
  stat.np$sortstate <- factor(stat.np$state, levels =  stat.np$state[order(-stat.np$num_cpg_total)])
  stat.np$Context <- "Non-promoter"
  stat.np$color <- "Red"
  
  color <- data.frame(color=c("Skyblue", "Red"), row.names=c("Promoter", "Non-promoter"))
  stat <- rbind(stat.p, stat.np)
  
  p <- ggplot(stat, aes(x = sortstate, y = num_cpg_total, fill=Context)) + 
    geom_bar(position="dodge") + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, colour='black'))
  p
  return(p)
}

plothist.stat_item <- function(stat, state, item, color, title, xlab, ylab,
                               cexaxis) {
  stat.select <- stat[state, ]
  stat.sorted <- stat.select[order(-stat.select[, item]), ]
  h <- max( c(stat.sorted[1, item] / 4, stat.sorted[2, item]) * 1.05 )
  p <- plot(1:length(state), stat.sorted[, item], ylim=c(0, h), col=color, ann=FALSE, xaxt='n', type="h")
  title(main=title)
  axis(1, at = 1:length(state), labels=rownames(stat.sorted), las=2, cex.axis=cexaxis)
  axis.break(axis=2, breakpos=h*0.9, pos=NA, bgcol="white", breakcol="black", style="slash",brw=0.02)
  mtext(xlab, side=1, line=3)  
  mtext(ylab, side=2, line=2)  
  return(p)
}



gridtab.samplelist <- function(list.sample, fs) {
  g <- grid.table(list.sample,
                  gpar.coretext = gpar(fontsize = fs, col="black"), 
                  gpar.coltext  = gpar(fontsize = fs,col="black", fontface="bold"), 
                  gpar.rowtext  = gpar(fontsize = fs, fontface="bold"))
  return(g)
}


plotstrip.stat_item <- function(stat, order, item, legendinsect, top=11, scale=1000, osep='', fontsize) {
  stat <- stat[order(-stat[, item]), ]
  span <- sum(stat[, item]) / scale
  span2 <- sum(stat[1:top, item]) / scale
  num_state <- dim(stat)[1]
  
  stat$color <- '#FFFFFF'
  stat.color <- stat[1:top, c('state', item, 'strip_pos')]
  stat.color$freq <- sapply(1:top, function(x) nnzero(strsplit(stat.color[x, 'state'], osep)[[1]]) )
  stat.color <- stat.color[order(-stat.color[, 'freq']), ]
  #colfunc <- colorRampPalette(c("blue", "white"))
  stat.color$color <- brewer.pal(top, "Spectral")
  stat[rownames(stat.color), 'color'] <- stat.color$color
  
  stat <- stat[order(stat[, 'strip_pos']), ]
  
  par(mar=c(5.1, 2, 4.1, 10), xpd=TRUE)
  plot(c(0, span), c(0, 6), type = "n", xlab="kb", ylab="", yaxt='n')
  text(span/2, 5.2, "All states", cex=1)
  cursor <- 0
  text(cursor, 5.2, rownames(stat[1, ]), cex=0.8)
  for (i in 1:num_state) {
    nextcursor <- cursor+stat[i, item]/scale
    rect(cursor, 4, nextcursor, 5, col=stat[i, 'color'], border="black")
    cursor <- nextcursor
  }
  text(cursor, 5.2, rownames(stat[num_state, ]), cex=0.8)
  
  stat.color <- stat.color[order(stat.color[, 'strip_pos']), ]
  text(span/2, 3.2, paste0("Top ", top, " states"), cex=1)
  cursor <- (span - span2) / 2
  text(cursor, 3.2, rownames(stat.color[1, ]), cex=0.8)
  for (i in 1:top) {
    nextcursor <- cursor+stat.color[i, item]/scale
    rect(cursor, 2, nextcursor, 3, col=stat.color[i, 'color'], border="black")
    cursor <- nextcursor
  }
  text(cursor, 3.2, rownames(stat.color[top, ]), cex=0.8)
  stat.color <- stat.color[order(-stat.color[, item]), ]
  legend("topright", inset=c(-legendinsect,0), legend=rownames(stat.color),
         fill=stat.color$color, cex=fontsize)
  stat.color <- stat.color[order(-stat.color[, 'freq']), ]
  text(400*span/scale, 1.2, rownames(stat.color[1, ]), cex=0.8)
  text(600*span/scale, 1.2, rownames(stat.color[top, ]), cex=0.8)
  color.legend(400*span/scale, 0.5, 600*span/scale, 1, legend='', rect.col=stat.color$color, gradient="x")
}

plotscatter.pos_item <- function(stat, posmat, item='num_cpg_total', top=20, dotsize=10) {
  vals <- unlist(t(posmat))
  vals <- vals[vals != 0]
  x <- data.frame(matrix(0, nrow <- length(vals), ncol <-3))
  colnames(x) <- c('state', 'pos', 'val')
  last <- 1
  for(state in rownames(posmat)) {
    v <- posmat[state, ]
    v <- v[state, v!=0]
    if (length(v) > 0) {
      end <- last + length(v) - 1
      x[last:end, 'state'] <- state
      x[last:end, 'pos'] <- as.numeric(colnames(v))
      x[last:end, 'val'] <- t(v[1, ])
      last <- end + 1 
    }
  }
  
  states.select <- rownames(stat[order(-stat[, item]), ])[1:top]
  x.select <- x[x[, 'state'] %in% states.select, ]
  
  p <- ggplot(x.select, aes(state, pos))
  p <- p + geom_point(aes(color = val), size = dotsize) + coord_flip() +
       scale_colour_gradient(low="lightskyblue1", high="brown")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
  return(p)
  #p + geom_point(aes(size = val)) + scale_size_continuous(range = c(1, maxplot))
}


multihist.tab_item <- function(stat, tab, item, sortitem, nrow=4, ncol=3, cutfirst=0.2) {
  stat <- stat[order(-stat[, sortitem]), ]
  state <- rownames(stat[stat[, sortitem] > 0, ])
  page_nfig <- nrow * ncol
  fig_left <- length(state)
  state_idx <- 1
  while(fig_left > 0) {
    par(mar=c(4, 4, 4, 4))
    layout(matrix(1:page_nfig, nrow = nrow, byrow = TRUE))
    for (i in state_idx:(state_idx+page_nfig-1)) {
      x <- tab[tab[, 'state']==state[i], item]
      largest <- max(x)
      v <- hist(x, breaks=100, freq=FALSE, plot=FALSE)
      y2 <- max(v$density[1] * cutfirst, v$density[2])
      hist(x, breaks=100, freq=FALSE, xlim=c(0, largest), ylim=c(0, y2),
           xlab=item, main=paste0(item, " of state ", state[i]))
      axis.break(axis=2, breakpos=0.8*y2, pos=NA, bgcol="white", breakcol="black", style="slash",brw=0.02)
    }
    state_idx <- state_idx+page_nfig-1
    fig_left <- fig_left - page_nfig
  }
}

#plot.stat_pie <- function()
#  pdf(file=file.path(dir, "piechart.pdf"), width=10, height=6, pointsize=10)
#layout(matrix(1::dim(list.sample)[1], 3, 2, byrow = TRUE))
#for(i in 1:dim(list.sample)[1]) {
#  state.select <- state[sapply(state, function(x) strsplit(x, "")[[1]][i] == 1)]
#  val <- data.frame(matrix(0, nrow=length(state.select), ncol=1))
#  colnames(val) <- 'total_length'
#  val$state <- state.select
#  rownames(val) <- state.select
#  for(s in state.select) {
#    val[s, 'total_length'] <- sum(tab[tab[,'state']==s, 'size'])
#  }
#  val <- val[order(-val[, 'total_length']), ]
#  pie(val[, 'total_length'], labels = NA, col = rainbow(length(state.select)),
#      main=paste0("Consider sample ", i))
#  legend("left", rownames(val), cex=0.5, fill=rainbow(length(state.select)) )
#}

