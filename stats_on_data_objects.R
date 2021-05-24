library(tidyverse)

tils <- readRDS('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/Seurat_Objects/')
tils@meta.data$Row.names <- rownames(tils@meta.data)
mta <- tils@meta.data
mta$Patient <- mta$patient
mta$Patient <- ifelse(mta$Patient=='p11' & grepl('rel', mta$sample), 'p11_rel', mta$Patient)
mta$Patient <- ifelse(mta$Patient=='p11' & grepl('pre', mta$sample), 'p11_pre', mta$Patient)
cat('number of cells per population:')
table(mta$Patient)

p2.rownames <- mta[mta$Patient == 'p2',]$Row.names
p11pre.rownames <- mta[mta$Patient == 'p11_pre',]$Row.names
p11rel.rownames <- mta[mta$Patient == 'p11_rel',]$Row.names
p6.rownames <- mta[mta$Patient == 'p6',]$Row.names
p15.rownames <- mta[mta$Patient == 'p15',]$Row.names


count.mean = mean(tils$nCount_RNA)
count.median = median(tils$nCount_RNA)
count.sd = sd(tils$nCount_RNA)

count.mean; count.median; count.sd


totcounts = colSums(tils@assays$RNA@counts)
countspergene = rowSums(tils@assays$RNA@counts!=0)
numgenes = colSums(tils@assays$RNA@counts!=0)

p2.totcounts <- totcounts[p2.rownames]; length(p2.totcounts)
p6.totcounts <- totcounts[p6.rownames]; length(p6.totcounts)
p15.totcounts <- totcounts[p15.rownames]; length(p15.totcounts)
p11pre.totcounts <- totcounts[p11pre.rownames]; length(p11pre.totcounts)
p11rel.totcounts <- totcounts[p11rel.rownames]; length(p11rel.totcounts)


p2.numgenes <- numgenes[p2.rownames]; length(p2.numgenes)
p6.numgenes <- numgenes[p6.rownames]; length(p6.numgenes)
p15.numgenes <- numgenes[p15.rownames]; length(p15.numgenes)
p11pre.numgenes <- numgenes[p11pre.rownames]; length(p11pre.numgenes)
p11rel.numgenes <- numgenes[p11rel.rownames]; length(p11rel.numgenes)


allstats <- getStats(p2.numgenes, 'p2_genecount')
allstats <- rbind(allstats, getStats(p2.totcounts, 'p2_total_counts_per_cell'))
allstats <- rbind(allstats, getStats(p15.numgenes, 'p15_genecount'))
allstats <- rbind(allstats, getStats(p15.totcounts, 'p15_total_counts_per_cell'))
allstats <- rbind(allstats, getStats(p6.numgenes, 'p6_genecount'))
allstats <- rbind(allstats, getStats(p6.totcounts, 'p6_total_counts_per_cell'))
allstats <- rbind(allstats, getStats(p11pre.numgenes, 'p11pre_genecount'))
allstats <- rbind(allstats, getStats(p11pre.totcounts, 'p11pre_total_counts_per_cell'))
allstats <- rbind(allstats, getStats(p11rel.numgenes, 'p11rel_genecount'))
allstats <- rbind(allstats, getStats(p11rel.totcounts, 'p11rel_total_counts_per_cell'))


getStats <- function(numbers, name) {
  mn <- mean(numbers)
  med <- median(numbers)
  std <- sd(numbers)
  min <- min(numbers)
  max <- max(numbers)
  
  row <- data.frame('mean'=mn, 'median'=med, 'std'=std, 'min'=min, 'max'=max)
  rownames(row) <- name
  return(row)
}

p2.cellstuff <- data.frame(matrix(ncol=length(p2.numgenes), nrow=2))
colnames(p2.cellstuff) <- p2.rownames; rownames(p2.cellstuff) <- c('gene_count', 'total_counts')
p2.cellstuff['gene_count',] <- p2.numgenes
p2.cellstuff['total_counts',] <- p2.totcounts

p6.cellstuff <- data.frame(matrix(ncol=length(p6.numgenes), nrow=2))
colnames(p6.cellstuff) <- p6.rownames; rownames(p6.cellstuff) <- c('gene_count', 'total_counts')
p6.cellstuff['gene_count',] <- p6.numgenes
p6.cellstuff['total_counts',] <- p6.totcounts

p15.cellstuff <- data.frame(matrix(ncol=length(p15.numgenes), nrow=2))
colnames(p15.cellstuff) <- p15.rownames; rownames(p15.cellstuff) <- c('gene_count', 'total_counts')
p15.cellstuff['gene_count',] <- p15.numgenes
p15.cellstuff['total_counts',] <- p15.totcounts

p11pre.cellstuff <- data.frame(matrix(ncol=length(p11pre.numgenes), nrow=2))
colnames(p11pre.cellstuff) <- p11pre.rownames; rownames(p11pre.cellstuff) <- c('gene_count', 'total_counts')
p11pre.cellstuff['gene_count',] <- p11pre.numgenes
p11pre.cellstuff['total_counts',] <- p11pre.totcounts

p11rel.cellstuff <- data.frame(matrix(ncol=length(p11rel.numgenes), nrow=2))
colnames(p11rel.cellstuff) <- p11rel.rownames; rownames(p11rel.cellstuff) <- c('gene_count', 'total_counts')
p11rel.cellstuff['gene_count',] <- p11rel.numgenes
p11rel.cellstuff['total_counts',] <- p11rel.totcounts


write.table(p2.cellstuff, '~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/cell_stats/patient2_cellstats.txt', sep='\t')
write.table(p6.cellstuff, '~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/cell_stats/patient6_cellstats.txt', sep='\t')
write.table(p15.cellstuff, '~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/cell_stats/patient15_cellstats.txt', sep='\t')
write.table(p11pre.cellstuff, '~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/cell_stats/patient11pre_cellstats.txt', sep='\t')
write.table(p11rel.cellstuff, '~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/cell_stats/patient11rel_cellstats.txt', sep='\t')
write.table(allstats, '~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/tils_combined/cell_stats/all_stats.csv', sep=',')

