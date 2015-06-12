###############################################################################
dir.data <- '~/data/subHMR/state_stat/0520/brain_blood_33_2'
core <- "111,222"
file.tab <- file.path(dir.data, 's200.bed')
file.idlist <- file.path(dir.data, 's200.list')
#file.category <- '~/data/hmr/all/category.txt'
file.category <- ''
file.samplelist <- '~/data/hmr/all/hmrID.txt'
###############################################################################
source("~/codes/HMRbase/rscripts/lib/subHMRfiles.R")
list.sample <- loadSampleList(file.samplelist, file.idlist, file.category)

tssfile <- paste0(file.tab, '.tss')
pfile <- paste0(file.tab, '.p')
tab <- loadSHMR_tss(tssfile, pfile, list.sample)
source("~/codes/HMRbase/rscripts/lib/state_proc.R")
idmap <- state.assign_tissue_id(tab$state, list.sample)
tab$state <- idmap[tab$state, 'state']

###############################################################################

tab.p <- tab[tab[, 'promoter'], ]
tab.np <- tab[!tab[, 'promoter'], ]

state.p <- unique(tab.p$state)
stat.p <- data.frame(count=matrix(numeric(), nrow=length(state.p)), row.names=state.p)
stat.p <- stat.count(tab.p, stat.p)
stat.p <- stat.sum(tab.p, stat.p, 'size')
stat.p <- stat.sum(tab.p, stat.p, 'num_cpg')

state.np <- unique(tab.np$state)
stat.np <- data.frame(count=matrix(numeric(), nrow=length(state.np)), row.names=state.np)
stat.np <- stat.count(tab.np, stat.np)
stat.np <- stat.sum(tab.np, stat.np, 'size')
stat.np <- stat.sum(tab.np, stat.np, 'num_cpg')

source("~/codes/HMRbase/rscripts/lib/state_plot.R")
pdf(file=paste0(dir.data, '/fig/sample.pdf'), width=10, height=10, pointsize=10)
par(mar=c(6,5,3,5))
gridtab.samplelist( as.data.frame( t(list.sample[, c('sample_id','project','tissue','tissue_id')]) ),
                    7)
dev.off()

pdf(file=paste0(dir.data, '/fig/size_total.pdf'), width=10, height=10, pointsize=10)
plothist.size.p(stat.p, stat.np)
dev.off()
pdf(file=paste0(dir.data, '/fig/cpg_total.pdf'), width=10, height=10, pointsize=10)
plothist.cpg.p(stat.p, stat.np)
dev.off()


pdf(file=paste0(dir.data, '/fig/size_total_20.pdf'), width=10, height=10, pointsize=10)
plothist.size.p(stat.p, stat.np)
dev.off()
pdf(file=paste0(dir.data, '/fig/cpg_total_20.pdf'), width=10, height=10, pointsize=10)
plothist.cpg.p(stat.p, stat.np)
dev.off()


source("~/codes/HMRbase/rscripts/lib/state_proc.R")
pos.p <- pos.stateorder(tab.p, stat.p, item='size', midstate=core)
posmat.p <- pos.p[[1]]
order.p <- pos.p[[2]]
pos.np <- pos.stateorder(tab.np, stat.np, item='size', midstate=core)
posmat.np <- pos.np[[1]]
order.np <- pos.np[[2]]

#par(mar=c(6,5,3,5))
#layout(matrix(c(1, 1, 2, 3), nrow = 4, byrow = TRUE))
#plothist.stat_item(stat, rownames(stat), 'num_cpg_total', 'black', "total # of CpG", NA, 'number', 1)
#plotstrip.stat_item(stat, order, 'num_cpg_total', 0.08, top=11, scale=1000, fontsize=0.9)
#pushViewport(viewport(x=.5, y=.15, height=.4))
#gridtab.samplelist( as.data.frame( t(list.sample[, c('sample_id','project','tissue','tissue_id')]) ),
#                    7)
#popViewport()
source("~/codes/HMRbase/rscripts/lib/state_plot.R")
pdf(file=paste0(dir.data, '/fig/order.pdf'), width=10, height=10, pointsize=10)
par(mar=c(6,5,3,5))
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
plotscatter.pos_item(stat.p, posmat.p, 'size_total')
plotscatter.pos_item(stat.np, posmat.np, 'size_total')
dev.off()
