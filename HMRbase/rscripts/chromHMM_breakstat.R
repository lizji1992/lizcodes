library(ggplot2)

detect_break <- function(tab) {
  for (i in 1:degree) {
    num <- dim(tab)[1]
    tab$breaks <- FALSE
    tab[2:(num-1), 'breaks'] <- apply(2:(num-1),
                                       function(x)
                                         tab[x, 'num_member'] > 1 &&
                                         tab[(x-1), 'state'] == tab[(x+1), 'state'] &&
                                         tab[x, 'size'] < min(tab[(x-1), 'size'], tab[(x+1), 'size'])
    )
    
    tab$breaked <- FALSE
    index_breaks <- tab[tab[, 'breaks'], 'index']
    index_breaked_left <- index_breaks - matrix(1, nrow=length(index_breaks), ncol=1)
    index_breaked_right <- index_breaks + matrix(1, nrow=length(index_breaks), ncol=1)
    tab[index_breaked_left, 'breaked'] <- TRUE
    tab[index_breaked_right, 'breaked'] <- TRUE
  }
  return(tab)
}

loadBED_cpg <- function(tabfile) {
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c('character',rep('numeric',2), 'character', 'numeric'),
                    col.names=c('chr', 'start', 'end', 'state', 'num_cpg'))
  tab$size <- (tab$end - tab$start + matrix(1, nrow=dim(tab)[1], ncol=1))
  return(tab)
}

loadBED_break <- function(tabfile) {
  tab <- read.table(tabfile, header=FALSE, sep='\t',
                    colClasses=c(rep('NULL', 4), 'character'),
                    col.names=c(rep('', 4), 'break_type'))
  return(tab)
}
###############################################################################
dir.data <- '~/data/subHMR/subhmr_chromHMM/subhmr_chromHMM_cpg_s25'
file.cpgtab <- file.path(dir.data, 'subhmr_25_segments_cpg.bed')
step <- 3

dir.out <- file.path(dir.data, 'fig')
file.breaktab <- file.path(dir.data,
                           paste0('subhmr_25_segments_', step, 'step-break.bed'))

###############################################################################
tab <- loadBED_cpg(file.cpgtab)
tab <- cbind(tab, loadBED_break(file.breaktab))

###############################################################################


tab$cpg_density <- (tab$num_cpg / tab$size)

ttest <- t.test(tab[tab[, 'break_type'] == "broken", 'cpg_density'],
                tab[tab[, 'break_type'] == "break", 'cpg_density'], alternative='greater')


pdf(file=file.path(dir.out, paste0(step, 'step_break-cpgd-stats.pdf')), width=10, height=10, pointsize=10)
tab.break <- tab[tab[, 'break_type'] != "na", ]
#tab.break <- na.omit(tab[, c('cpg_density', 'break_type')])
p <- ggplot(data=tab.break, aes(tab.break$cpg_density, y=..count.., fill=tab.break$break_type)) +
     geom_density(alpha=0.3) + theme_bw()
p
dev.off()

sink(file.path(dir.out, paste0(step, 'step_break-cpgd-ttest.txt')))
ttest
sink()
#geom_density(alpha=0.4) + theme_bw() + coord_cartesian(xlim = c(-0.0001, 0.01))
#p <- p + annotate("text", x = 0.005, y = 2e6,
#                  label = paste0("Mean-breaked:", mean(tab[tab[, 'breaked'], 'cpg_density']), "\n",
#                                 "Mean-breaks:", mean(tab[tab[, 'breaks'], 'cpg_density']), "\n",
#                                 "p-val:", ttest$p.val, " breaks don't have significantly higher CpG density" ))


