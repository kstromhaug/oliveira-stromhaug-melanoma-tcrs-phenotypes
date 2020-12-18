library(SingleR)
library(Seurat)
library(scater)

wu = readRDS('wu_merged_cd8_withscores.rds')
tils = readRDS('tils_combinde_cd8_withscores.rds')
yost = readRDS('yost_justcd8_20200826.rds')

Idents(wu)<-wu$tcell.type
wu.nonaive <- subset(wu, idents='8.6-KLRB1', invert=T)
# markers = readRDS('wu_cluster_markers_list_sixrem.rds')

tils$final_cluster <- tils$seurat_clusters
tils$final_cluster <- as.numeric(as.character(tils$final_cluster))
tils$final_cluster <- ifelse(tils$final_cluster %in% c(0, 8, 11), 'term_ex', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 5, 'mitotic', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 4, 'prog_ex', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster %in% c(1,2), 'memory', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 10, 'NK_like_T', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 7, 'GD_like_T', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 3, 'CD8_eff', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 6, 'CD8_MT', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 9, 'CD8_T_reg', tils$final_cluster)
tils$final_cluster <- ifelse(tils$final_cluster == 12, 'CD8_Naive', tils$final_cluster)


genesintersect = intersect(rownames(tils), rownames(wu)); length(genesintersect)
genesintersect = intersect(genesintersect, rownames(yost)); length(genesintersect)

trimmed_mean <- function(row) {
	row = row[order(row)]
	ind = length(row) * 0.1
	mn = mean(row[(ind+1):(length(row)-ind)])
	return(mn)
}

prep_data_for_training <- function(counts, obj, cell_categories, meta_othercol) {
	Idents(obj)<-obj@meta.data[,cell_categories]

	tcell.types = unique(obj@meta.data[,cell_categories])

	sigs = data.frame(matrix(nrow=nrow(counts), ncol=1))
	rownames(sigs) <- rownames(counts)
	colnames(sigs) <- 'gene'
	sigs$gene <- rownames(sigs)

	show('making SummarizedExperiment object')
	tcolD = obj@meta.data[,c(cell_categories, meta_othercol)]
	sumobj <- SummarizedExperiment(list(counts=counts), colData=tcolD)

	sumobj <- scater::logNormCounts(sumobj)
	logcounts <- sumobj@assays@data$logcounts
	show('normalized data, making lists')
	## interate through the cell types
	for (cell in tcell.types) {
		print(cell)

		expsub <- subset(obj, idents=cell)
		logcountssub <- logcounts[,colnames(expsub)]

		show('applying trimmed means')
		trms <- data.frame('trms'=apply(logcountssub, 1, trimmed_mean))
		colnames(trms)<-c(cell)

		sigs <- cbind(sigs, trms)
	}
	return(sigs)
}


prep_data_to_score <- function(obj, cell_categories, meta_othercol, genes_to_use) {
	objtest = as.matrix(obj@assays$RNA@counts[rownames(obj@assays$RNA@counts) %in% genes_to_use,])
	dim(objtest)
	tcolD = obj@meta.data[,c(cell_categories, meta_othercol)]
	objtest <- SummarizedExperiment(list(counts=objtest), colData=tcolD)
	objtest <- scater::logNormCounts(objtest)

	return(objtest)
}

wuexp = wu@assays$RNA@counts
wuexp = wuexp[genesintersect,]
drop = which(apply(wuexp, 1, max)<1)
wuexp = wuexp[-drop,]; dim(wuexp) ## 15887 genes left now
# wuexp = wuexp[intersect(rownames(wuexp),genesintersect),]

wuexp.n = wu.nonaive@assays$RNA@counts
wuexp.n = wuexp.n[genesintersect,]
drop = which(apply(wuexp.n, 1, max)<1)
wuexp.n = wuexp.n[-drop,]; dim(wuexp.n)

tilsexp = tils@assays$RNA@counts
tilsexp = tilsexp[genesintersect,]
drop = which(apply(tilsexp, 1, max)<1)
tilsexp = tilsexp[-drop,]; dim(tilsexp) ## 15225 genes left now

yostexp = yost@assays$RNA@counts
yostexp = yostexp[genesintersect,]
drop = which(apply(yostexp, 1, max)<1)
yostexp = yostexp[-drop,]; dim(yostexp) ## 13608 genes left now

finalgenes = intersect(rownames(tilsexp), rownames(wuexp)); length(finalgenes)
finalgenes = intersect(finalgenes, rownames(yostexp)); length(finalgenes)
## 13365 genes total to use now

wuexp = wuexp[finalgenes, ]
wuexp.n = wuexp.n[finalgenes, ]
tilsexp = tilsexp[finalgenes, ]
yostexp = yostexp[finalgenes, ]

### for each cluster in wu, find the 'robust centroids' that wu says they did for yost and the other dataset.



wusigs <- prep_data_for_training(wuexp, wu, 'tcell.type', 'meta.barcode')
write.table(wusigs, 'wusigs_training2_20200827.txt', sep='\t')
wusigs$gene <- NULL
tcell.types <- colnames(wusigs)
wutrained = trainSingleR(ref=as.matrix(wusigs), labels=tcell.types)

### train on Yost
yostsigs <- prep_data_for_training(yostexp, yost, 'cluster', 'Row.names')
write.table(yostsigs, 'yostsigs_training2_20200827.txt', sep='\t')
yostsigs$gene <- NULL
tcell.types <- colnames(yostsigs)
yosttrained = trainSingleR(ref=as.matrix(yostsigs), labels=tcell.types)

###################### SCORE THE DATASETS ########################

## LABEL TILS WITH WU LABELS
scdtest <- prep_data_to_score(tils, 'final_cluster', 'Row.names', finalgenes)
tils.wuscored <- classifySingleR(test=scdtest, trained=wutrained, fine.tune=TRUE)
colnames(tils.wuscored) <- paste0('wutrained.', colnames(tils.wuscored))

## LABEL TILS WITH YOST LABELS
tils.yostscored <- classifySingleR(test=scdtest, trained=yosttrained, fine.tune=TRUE)
colnames(tils.yostscored) <- paste0('yosttrained.', colnames(tils.yostscored))

##################### SAVE NECESSARY DATA FILES

tils.labels <- cbind(tils.wuscored[,c('wutrained.labels', 'wutrained.pruned.labels')], 
					tils.yostscored[,c('yosttrained.labels', 'yosttrained.pruned.labels')])

write.table(tils.labels, 'tils_with_wuyost_labels_finetuned_20200827.txt', sep='\t')

saveRDS(tils, 'tils_singlerd.rds')




