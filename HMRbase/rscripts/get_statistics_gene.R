source("~/codes/HMRbase/rscripts/lib/subHMRfiles.R")
source("~/codes/HMRbase/rscripts/lib/state_proc.R")
source("~/codes/HMRbase/rscripts/lib/state_plot.R")
library(car) 

###############################################################################
dir.data <- '~/data/subHMR/state_stat/0427/1'
file.tab <- file.path(dir.data, 'all.bed')
file.idlist <- file.path(dir.data, 'all.list')

dir.gene <- '~/data/subHMR/examples/updated_dataset'
file.genelist <- file.path(dir.gene, 'genetable.txt')

file.samplelist <- '~/data/hmr/all/hmrID.txt'
###############################################################################

tab <- loadSHMR(file.tab)
list.sample <- loadSampleList(file.samplelist, file.idlist, file.category)
list.gene <- loadIntervals_UCSC(file.genelist)
#------------------------------------------------------------------------------
pdf(file=file.path(dir.data, 'gene_states.pdf'), width=10, height=10, pointsize=10)
par(mar=c(6,5,3,5))
gridtab.samplelist( list.sample[, c('sample_id','project','tissue','tissue_id')], fs=7)
for (gene in rownames(list.gene)) {
  tab.select <- hmrs.in_range(tab, list.gene[gene, ])

  state <- unique(tab.select$state)
  stat <- data.frame(count=matrix(numeric(), nrow=length(state)), row.names=state)
  stat$state <- state
  stat <- stat.size(tab.select, stat)
  rownames(stat) <- state.assign_tissue_id(stat$state, list.sample, osep=',', tsep=',')
  stat$strip_pos <- strip_pos.state(stat$state)
  
  par(mar=c(6,5,3,5))
  layout(matrix(c(1,1,2,3), nr = 4, byrow = TRUE))
  plothist.stat_item(stat, rownames(stat), 'len_total', 'black',
                     paste0("total length gene: ", gene), NA, 'length(bp)', 1)
  plot.new()
  plotstrip.stat_item(stat, 'len_total', 0.1, top=11, scale=1000, fontsize=0.6)
}
dev.off()


#plot.stat(stat, file.path(dir.data, "state_stat.pdf"), state, list.sample, state)  
#plot(tab[, "state"], tab[, "size"], main="Scatterplot Example", 
#     xlab="State", ylab="Size(bp)", pch=19)
#scatterplot(state, tab[, "size"])
#state[sapply(state, function(x) count(strsplit(x, "")[[1]])[2,2] > 1)]



