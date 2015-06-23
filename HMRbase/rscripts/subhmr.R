library(ggplot2)

map.state_id <- function(file.id) {
 id <- read.table(file.id, header=FALSE, sep='\t',
                  colClasses=c('character', 'character'),
                  col.names=c('id', 'state'))
  rownames(id) <- id$id
  return(id)
}

###############################################################################

dir.data <- "~/data/subHMR/subhmr/chromHMM/choose_ct_cpg_s25"
pref <- file.path(dir.data, 'bed/subhmr_25_segments_cpg.bed')
file.emi <- file.path(dir.data, 'emissions_25.txt')

source("~/codes/HMRbase/rscripts/lib/BedFiles.R")
id <- map.state_id(paste0(pref, '.id')) 
p <- loadBED_p(paste0(pref, '.p'))
np <- loadBED_np(paste0(pref, '.np'))

#################################################

p$loci <- 'p'
np$loci <- 'np'
tab.p_np <- rbind(p[, c('chr', 'start', 'end', 'id', 'num_cpg', 'loci')],
                  np[, c('chr', 'start', 'end', 'id', 'num_cpg', 'loci')])
rm(p)
rm(np)
# count size
tab.p_np$size <- (tab.p_np$end - tab.p_np$start)
tab.p_np.g200 <- tab.p_np[tab.p_np[, 'size']>200, ]

pdf(file=file.path(dir.data, 'fig', 'sizes_pvsnp.pdf'), width=10, height=6, pointsize=10)
p <- ggplot(data=tab.p_np.g200, aes(tab.p_np.g200$size, y=..count.., fill=tab.p_np.g200$loci)) +
     geom_density(alpha=0.3) + theme_bw()
p
dev.off()

#################################################
source("~/codes/HMRbase/rscripts/lib/chromHMM_utils.R")
emi_matrix <- emission.load(file.emi)
emi <- cbind(emi_matrix, sum_emi=rowSums(emi_matrix))
tab.p_np$sum_emi <- emi[id[tab.p_np$id, 'state'], 'sum_emi']

source("~/codes/HMRbase/rscripts/lib/plot_utils.R")
pdf(file=file.path(dir.data, 'fig', "sum-emission_pvsnp.pdf"),width=10, height=6, pointsize=10)
tab.p <- tab.p_np[tab.p_np[, 'loci']=='p', ]
p1 <- ggplot(data=tab.p, aes(tab.p$sum_emi, y=..count..)) + geom_density(fill='skyblue') + theme_bw()
tab.np <- tab.p_np[tab.p_np[, 'loci']=='np', ]
p2 <- ggplot(data=tab.np, aes(tab.np$sum_emi, y=..count..)) + geom_density(fill='red') + theme_bw()
multiplot(p1, p2, cols=2)
dev.off()

#################################################
state.core <- emission.find_core(emi_matrix)
tab.p_np$state <- id[tab.p_np$id, 'state']
tab.core <- tab.p_np[tab.p_np[, 'state']==state.core, ]
pdf(file=file.path(dir.data, 'fig', 'core-size_pvsnp.pdf'), width=10, height=6, pointsize=10)
p <- ggplot(data=tab.core, aes(tab.core$size, y=..count.., fill=tab.core$loci)) +
  geom_density(alpha=0.3) + theme_bw()
p
dev.off()

#################################################



