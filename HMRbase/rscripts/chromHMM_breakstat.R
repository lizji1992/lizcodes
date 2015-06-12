library(ggplot2)

###############################################################################
dir.data <- '~/data/subHMR/subhmr_chromHMM/subhmr_chromHMM_cpg_s25'
file.tab <- file.path(dir.data, 'subhmr_25_segments_countcpg.bed')
#dir.data <- '~/data/subHMR/subhmr_chromHMM/b50_s25'
#file.tab <- file.path(dir.data, 'subhmr1_25_segments_countcpg.bed')

###############################################################################
source("~/codes/HMRbase/rscripts/lib/subHMRfiles.R")
tab <- loadSHMR_chromHMM(file.tab)

###############################################################################

source("~/codes/HMRbase/rscripts/lib/subHMR_analysis.R")
tab <- assign.group(tab)
tab$index <- as.numeric(rownames(tab))

# find breaks
num <- dim(tab)[1]
tab$breaks <- FALSE
tab[2:(num-1), 'breaks'] <- sapply(2:(num-1),
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


tab$cpg_density <- (tab$num_cpg / tab$size)
tab[tab[, 'breaks'], 'break_type'] <- 'breaks'
tab[tab[, 'breaked'], 'break_type'] <- 'breaked'
ttest <- t.test(tab[tab[, 'breaked'], 'cpg_density'], tab[tab[, 'breaks'], 'cpg_density'], alternative='greater')


pdf(file=paste0(dir.data, '/break_cpgdensity.pdf'), width=10, height=10, pointsize=10)
tab.break <- na.omit(tab[, c('cpg_density', 'break_type')])
p <- ggplot(data=tab.break, aes(tab.break$cpg_density, y=..count.., fill=tab.break$break_type)) +
     geom_density(alpha=0.4) + theme_bw()
p
dev.off()

sink(paste0(dir.data, '/break_ttest.txt'))
ttest
sink()
#geom_density(alpha=0.4) + theme_bw() + coord_cartesian(xlim = c(-0.0001, 0.01))
#p <- p + annotate("text", x = 0.005, y = 2e6,
#                  label = paste0("Mean-breaked:", mean(tab[tab[, 'breaked'], 'cpg_density']), "\n",
#                                 "Mean-breaks:", mean(tab[tab[, 'breaks'], 'cpg_density']), "\n",
#                                 "p-val:", ttest$p.val, " breaks don't have significantly higher CpG density" ))


