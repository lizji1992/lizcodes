###############################################################################
dir.data <- '~/data/subHMR/state_stat/0520/blood_embryo_33_1'
core <- "111,222"
file.tab <- file.path(dir.data, 's200.bed')
file.idlist <- file.path(dir.data, 's200.list')
#file.category <- '~/data/hmr/all/category.txt'
file.category <- ''
file.samplelist <- '~/data/hmr/all/hmrID.txt'
###############################################################################

source("~/codes/HMRbase/rscripts/lib/subHMRfiles.R")
list.sample <- loadSampleList(file.samplelist, file.idlist, file.category)
tab <- loadSHMR2(file.tab, list.sample)

###############################################################################

source("~/codes/HMRbase/rscripts/lib/state_proc.R")
state <- unique(tab$state)
stat <- data.frame(count=matrix(numeric(), nrow=length(state)), row.names=state)
stat <- stat.count(tab, stat)
stat <- stat.sum(tab, stat, 'size')
stat <- stat.sum(tab, stat, 'num_cpg')
stat <- stat.sum(tab, stat, 'score')
#tab$dist_left <- (tab$start - tab$istart)
#stat <- stat.mean(tab, stat, 'dist_left')
idmap <- state.assign_tissue_id(stat$state, list.sample)
stat$state <- idmap[stat$state, 'state']
rownames(stat) <- stat$state
tab$state <- idmap[tab$state, 'state']

source("~/codes/HMRbase/rscripts/lib/state_proc.R")
pos <- pos.stateorder(tab, stat, item='num_cpg', midstate=core)
posmat <- pos[[1]]
order <- pos[[2]]

source("~/codes/HMRbase/rscripts/lib/state_plot.R")
pdf(file=paste0(dir.data, '/fig/sample.pdf'), width=10, height=10, pointsize=10)
par(mar=c(6,5,3,5))
gridtab.samplelist( as.data.frame( t(list.sample[, c('sample_id','project','tissue','tissue_id')]) ),
                    7)
dev.off()

pdf(file=paste0(dir.data, '/fig/size_total.pdf'), width=10, height=10, pointsize=10)
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
plothist.stat_item(stat, rownames(stat), 'size_total', 'black', "total length", NA, 'length(bp)', 1)
pushViewport(viewport(x=.5, y=.15,height=.4))
gridtab.samplelist( as.data.frame( t(list.sample[, c('sample_id','project','tissue','tissue_id')]) ),
                    7)
dev.off()
pdf(file=paste0(dir.data, '/fig/cpg_total.pdf'), width=10, height=10, pointsize=10)
par(mar=c(6,5,3,5))
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
plothist.stat_item(stat, rownames(stat), 'num_cpg_total', 'black', "total CpGs", NA, 'length(bp)', 1)
pushViewport(viewport(x=.5, y=.15,height=.4))
gridtab.samplelist( as.data.frame( t(list.sample[, c('sample_id','project','tissue','tissue_id')]) ),
                    7)
dev.off()
pdf(file=paste0(dir.data, '/fig/score_total.pdf'), width=10, height=10, pointsize=10)
par(mar=c(6,5,3,5))
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
plothist.stat_item(stat, rownames(stat), 'score_total', 'black', "total score", NA, 'length(bp)', 1)
pushViewport(viewport(x=.5, y=.15,height=.4))
gridtab.samplelist( as.data.frame( t(list.sample[, c('sample_id','project','tissue','tissue_id')]) ),
                    7)
dev.off()

#par(mar=c(6,5,3,5))
#layout(matrix(c(1, 1, 2, 3), nrow = 4, byrow = TRUE))
#plothist.stat_item(stat, rownames(stat), 'num_cpg_total', 'black', "total # of CpG", NA, 'number', 1)
#plotstrip.stat_item(stat, order, 'num_cpg_total', 0.08, top=11, scale=1000, fontsize=0.9)
#pushViewport(viewport(x=.5, y=.15, height=.4))
#gridtab.samplelist( as.data.frame( t(list.sample[, c('sample_id','project','tissue','tissue_id')]) ),
#                    7)
#popViewport()
pdf(file=paste0(dir.data, '/fig/order.pdf'), width=10, height=10, pointsize=10)
par(mar=c(6,5,3,5))
layout(matrix(c(1, 2, 3), nrow = 3, byrow = TRUE))
plotscatter.pos_item(stat, posmat, 'size_total')
plotscatter.pos_item(stat, posmat, 'num_cpg_total')
plotscatter.pos_item(stat, posmat, 'score_total')
dev.off()
