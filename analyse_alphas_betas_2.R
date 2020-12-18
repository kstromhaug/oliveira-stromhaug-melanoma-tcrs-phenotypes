library(dplyr)
library(tidyverse)
library(ggpubr)
# library(Seurat)
library(openxlsx)
# library(pheatmap)

setwd("~/Dropbox (Partners HealthCare)/Melanoma p2 Analysis/")


#########################################################################################################
########################################## FOR THE PBMCS ################################################
#########################################################################################################

pbmc.cd8 = readRDS('pbmc_combined/pbmc.CD8.harmonized.20200418.rds')
pbmc.cd8 = readRDS('../pbmc_combined/pbmc.harmonized.20200415.rds')
pbmc.cd8$clusters <- pbmc.cd8$seurat_clusters

# assert(FALSE)
meta = pbmc.cd8@meta.data

cite.cols = colnames(meta)[grep('CITEseq.norm', colnames(meta))]; length(cite.cols)
rest.keep = c('patient', 'til.clonotype.family', 'TCR.Clone', 'sample', 'final.clonotype.family', 'final.family.counts', 'clusters', 'cd8.adt', 'cd4.adt', 'dp.adt', 'dn.adt', cite.cols)

p2.pbmc = meta %>% subset(patient == 'p2'); dim(p2.pbmc)
setdiff(rest.keep, colnames(p2.pbmc))
p2.pbmc = p2.pbmc[,rest.keep]
p2.origs = unique(p2.pbmc$sample); p2.origs
p2.labels = c('p2.post', 'p2.pre', 'p2.relapse')
p2.pbmc$patient.t <- 0
for (l in 1:length(p2.origs)) {
  p2.pbmc[p2.pbmc$sample==p2.origs[l],]$patient.t <- p2.labels[l]
}

p6.pbmc = meta %>% subset(patient=='p6' | patient == 'p11'); dim(p6.pbmc)
p6.pbmc = p6.pbmc[,rest.keep]
p6.origs = unique(p6.pbmc$sample); p6.origs
p6.pbmc$patient.t <- p6.pbmc$patient

p15.pbmc = meta %>% subset(patient=='p15'); dim(p15.pbmc)
p15.pbmc = p15.pbmc[,rest.keep]
p15.origs = unique(p15.pbmc$sample); p15.origs
p15.labels = c('p15.post', 'p15.pre')
p15.pbmc$patient.t <- 0
for (l in 1:length(p15.origs)) {
  p15.pbmc[p15.pbmc$sample==p15.origs[l],]$patient.t <- p15.labels[l]
}


pbmc.allmeta = rbind(p2.pbmc, p6.pbmc, p15.pbmc); dim(pbmc.allmeta)
pbmc.allmeta$cell.barcode <- rownames(pbmc.allmeta)

unique(pbmc.allmeta$patient)
unique(pbmc.allmeta$patient.t)

write.table(pbmc.allmeta, 'pbmc_combined/pbmc_cd8_metadata_0504.txt', sep='\t', row.names=F)


pbmc.allmeta = read.delim('pbmc_combined/pbmc_all_metadata_0504.txt', sep='\t')
cite.cols = colnames(pbmc.allmeta)[grep('CITEseq.norm', colnames(pbmc.allmeta))]; length(cite.cols)
#########################################################################################################
####################################### FOR THE TILS.CD8 ################################################
#########################################################################################################

tils.cd8 = readRDS('tils_combined/Seurat_Objects/tils.CD4.harmonized.20200514.rds')
tils.cd8 = readRDS('tils_combined/Seurat_Objects/tils.CD8.R0.6.harmonized.20200422.rds')
tils.cd8$clusters = tils.cd8$seurat_clusters
meta = tils.cd8@meta.data; dim(meta)
## SECELCT COLUMNS THAT WE WANT TO KEEP
cite.cols = colnames(meta)[grep('CITEseq.norm', colnames(meta))]; length(cite.cols)
rest.keep = c('patient', 'til.clonotype.family', 'TCR.Clone', 'sample', 'final.clonotype.family', 'final.family.counts', 'clusters', 'cd8.adt', 'cd4.adt', 'dp.adt', 'dn.adt', cite.cols)

## ANNOTATE THE PATIENT CLONOTYPES BY ADDING PATIENT IDENTITY TO CLONOTYPE NUMBER
firsthree = meta %>% subset(patient %in% c('p2', 'p6', 'p15')); dim(firsthree)
firsthree = firsthree[,rest.keep]

## TREAT SPECIAL CASE OF P11, WHICH HAS TWO POPULATION OF TILS
p11_pre = meta %>% subset(patient == 'p11' & sample %in% c('TIL_pre_CD45pos_CD3pos-1','TIL_pre_CD45pos_CD3pos-2','TIL_pre_CD45pos_CD3pos-3', 'TIL_CD45pos_CD3neg')); dim(p11_pre)
p11_pre$til.clonotype.family<-NULL
p11_pre$patient <- paste0(p11_pre$patient, '.pre')
p11_pre = p11_pre[,setdiff(colnames(p11_pre), 'til.rel.clonotype.family')]
names(p11_pre)[names(p11_pre)=='til.pre.clonotype.family']<-'til.clonotype.family'
rownames(p11_pre) <- gsub('p11', 'p11.pre', rownames(p11_pre))
p11_pre = p11_pre[,rest.keep]
p11_rel = meta %>% subset(patient == 'p11' & sample %in% c('TIL_rel-1', 'TIL_rel-2')); dim(p11_rel)
p11_rel$til.clonotype.family<-NULL
p11_rel$patient <- paste0(p11_rel$patient, '.rel')
p11_rel = p11_rel[,setdiff(colnames(p11_rel), 'til.pre.clonotype.family')]
names(p11_rel)[names(p11_rel)=='til.rel.clonotype.family']<-'til.clonotype.family'
rownames(p11_rel) <- gsub('p11', 'p11.rel', rownames(p11_rel))
p11_rel = p11_rel[,rest.keep]

## COMBINE THE THREE PATIENTS AGAIN, RENAME CLONOTYPES BY PATIENT
tils.allmeta = rbind(firsthree, p11_pre, p11_rel); dim(tils.allmeta); dim(meta)
tils.allmeta$patient.clone.meta <- paste0(tils.allmeta$patient, '-', tils.allmeta$til.clonotype.family)
tils.allmeta$cell.barcode <- rownames(tils.allmeta)
just_p11 = meta %>% subset(patient=='p11')

write.table(tils.allmeta, 'tils_combined/tils_cd4_updated_metadata_0515.txt', sep='\t')
cat('finished loading datatypes')
unique.clones.meta = unique(tils.allmeta$patient.clone.meta); length(unique.clones.meta)

# tils.allmeta = read.delim( '../tils_combined/tils_updated_metadata_0422.txt', sep='\t')

#########################################################################################################
######################################## FOR THE TILS ALL ###############################################
#########################################################################################################

tils = readRDS('tils_combined/Seurat_Objects/tils.harmonized.20200403.rds')
tils$clusters = tils$seurat_clusters
meta = tils@meta.data; dim(meta)
## SECELCT COLUMNS THAT WE WANT TO KEEP
cite.cols = colnames(meta)[grep('CITEseq.norm', colnames(meta))]; length(cite.cols)
rest.keep = c('patient', 'til.clonotype.family', 'TCR.Clone', 'sample', 'final.clonotype.family', 'final.family.counts', 'clusters', 'cd8.adt', 'cd4.adt', 'dp.adt', 'dn.adt', cite.cols)

## ANNOTATE THE PATIENT CLONOTYPES BY ADDING PATIENT IDENTITY TO CLONOTYPE NUMBER
firsthree = meta %>% subset(patient %in% c('p2', 'p6', 'p15')); dim(firsthree)
firsthree = firsthree[,rest.keep]

## TREAT SPECIAL CASE OF P11, WHICH HAS TWO POPULATION OF TILS
p11_pre = meta %>% subset(patient == 'p11' & sample %in% c('TIL_pre_CD45pos_CD3pos-1','TIL_pre_CD45pos_CD3pos-2','TIL_pre_CD45pos_CD3pos-3', 'TIL_CD45pos_CD3neg')); dim(p11_pre)
p11_pre$til.clonotype.family<-NULL
p11_pre$patient <- paste0(p11_pre$patient, '.pre')
p11_pre = p11_pre[,setdiff(colnames(p11_pre), 'til.rel.clonotype.family')]
names(p11_pre)[names(p11_pre)=='til.pre.clonotype.family']<-'til.clonotype.family'
rownames(p11_pre) <- gsub('p11', 'p11.pre', rownames(p11_pre))
p11_pre = p11_pre[,rest.keep]
p11_rel = meta %>% subset(patient == 'p11' & sample %in% c('TIL_rel-1', 'TIL_rel-2')); dim(p11_rel)
p11_rel$til.clonotype.family<-NULL
p11_rel$patient <- paste0(p11_rel$patient, '.rel')
p11_rel = p11_rel[,setdiff(colnames(p11_rel), 'til.pre.clonotype.family')]
names(p11_rel)[names(p11_rel)=='til.rel.clonotype.family']<-'til.clonotype.family'
rownames(p11_rel) <- gsub('p11', 'p11.rel', rownames(p11_rel))
p11_rel = p11_rel[,rest.keep]

## COMBINE THE THREE PATIENTS AGAIN, RENAME CLONOTYPES BY PATIENT
tils.allmeta = rbind(firsthree, p11_pre, p11_rel); dim(tils.allmeta); dim(meta)
tils.allmeta$patient.clone.meta <- paste0(tils.allmeta$patient, '-', tils.allmeta$til.clonotype.family)
tils.allmeta$cell.barcode <- rownames(tils.allmeta)
just_p11 = meta %>% subset(patient=='p11')

write.table(tils.allmeta, 'tils_combined/tils_all_updated_metadata_0429.txt', sep='\t')
cat('finished loading datatypes')
unique.clones.meta = unique(tils.allmeta$patient.clone.meta); length(unique.clones.meta)

tils.allmeta = read.delim( 'tils_combined/tils_updated_metadata_0422.txt', sep='\t')

#########################################################################################################
#################################### PREPARE CLONOTYPE FAMILY DOC #######################################
#########################################################################################################


### pull together all the clonotype match documents
p2.groups = read.xlsx('patient2/matches_tils_pbmcs_7_revised.xlsx')
p6.groups = read.xlsx('patient6/R_output/p6_clonotype_matches_revised.xlsx')
p15.groups = read.xlsx('patient15/p15_matches_tils_pbmcs-revised.xlsx')
p11.groups = read.xlsx('patient11/R_output/p11_clonotype_matches_revised.xlsx')
p11.groups.1 = p11.groups[p11.groups$origin %in% c('til_pre_CD45pos_CD3pos_1','til_pre_CD45pos_CD3pos_2','til_pre_CD45pos_CD3pos_3', 'til_CD45pos_CD3neg'), setdiff(colnames(p11.groups), c('til.rel.clonotype.family', 'til.post.counts'))]; dim(p11.groups.1)
names(p11.groups.1)[names(p11.groups.1)=='til.pre.clonotype.family']<-'til.clonotype.family'
names(p11.groups.1)[names(p11.groups.1)=='til.pre.counts']<-'til.counts'
p11.groups.2 = p11.groups[p11.groups$origin %in% c('til_rel_1', 'til_rel_2'), setdiff(colnames(p11.groups), c('til.pre.clonotype.family', 'til.pre.counts'))]; dim(p11.groups.2)
names(p11.groups.2)[names(p11.groups.2)=='til.rel.clonotype.family']<-'til.clonotype.family'
names(p11.groups.2)[names(p11.groups.2)=='til.post.counts']<-'til.counts'

nrow(p11.groups); nrow(p11.groups.1)+nrow(p11.groups.2)
# p11.groups.1.cells = p11.groups.1 %>% subset(is.na(til.clonotype.family)); dim(p11.groups.1.na)

keep = intersect(colnames(p2.groups), intersect(colnames(p6.groups), intersect(colnames(p11.groups.1), colnames(p15.groups)))); keep

p2.groups = p2.groups[,keep]; p2.groups$patient='p2'; p2.groups$cell.barcode<-paste0('p2_', p2.groups$cell.barcode)
p6.groups = p6.groups[,keep]; p6.groups$patient='p6'; p6.groups$cell.barcode<-paste0('p6_', p6.groups$cell.barcode)
p11.groups.1 = p11.groups.1[,keep]; p11.groups.1$patient='p11'; p11.groups.1$cell.barcode<-paste0('p11.pre_', p11.groups.1$cell.barcode)
p11.groups.1$patient <- paste0(p11.groups.1$patient, '.pre')
p11.groups.2 = p11.groups.2[,keep]; p11.groups.2$patient='p11'; p11.groups.2$cell.barcode<-paste0('p11.rel_', p11.groups.2$cell.barcode)
p11.groups.2$patient <- paste0(p11.groups.2$patient, '.rel')
p15.groups = p15.groups[,keep]; p15.groups$patient='p15'; p15.groups$cell.barcode<-paste0('p15_', p15.groups$cell.barcode)
all.families = rbind(p2.groups, p6.groups, p11.groups.1, p11.groups.2, p15.groups)
dim(all.families)
all.families$tcell.clone <- paste0(all.families$patient, all.families$clonotype); dim(all.families)
all.families$patient.clone <- paste0(all.families$patient, '-', all.families$til.clonotype.family)

write.table(all.families, 'all_clonotype_matches_20200424.txt', sep='\t')
all.families = read.delim('all_clonotype_matches_20200424.txt', sep='\t')
unique.clones = unique(all.families$patient.clone); length(unique.clones)
unique.clones = setdiff(unique.clones, c('p2-NA', 'p6-NA', 'p11.pre-NA', 'p11.rel-NA', 'p15-NA')); length(unique.clones)
length(unique.clones)
dim(all.families)

all.families <- AddCloneCategory(all.families)
not_found = all.families[!(all.families$cell.barcode %in% with.categ$cell.barcode),]; dim(not_found)

#########################################################################################################
############################## PREPARE CLONOTYPE FAMILY DOC FOR PBMCS ###################################
#########################################################################################################
p2.groups = read.xlsx('patient2/matches_tils_pbmcs_7_revised.xlsx')
p6.groups = read.xlsx('patient6/R_output/p6_clonotype_matches_revised.xlsx')
p15.groups = read.xlsx('patient15/p15_matches_tils_pbmcs-revised.xlsx')
p11.groups = read.xlsx('patient11/R_output/p11_clonotype_matches_revised.xlsx')
p2.groups$patient <- 'p2'
p6.groups$patient <- 'p6'
p11.groups$patient <- 'p11'
p15.groups$patient <- 'p15'

keep = intersect(colnames(p2.groups), intersect(colnames(p6.groups), intersect(colnames(p11.groups), colnames(p15.groups)))); keep
all.families.pbmc = rbind(p2.groups[,keep], p6.groups[,keep], p11.groups[,keep], p15.groups[,keep]); dim(all.families.pbmc)
all.families.pbmc$patient.t <- all.families.pbmc$patient
all.families.pbmc[all.families.pbmc$patient=='p2' & all.families.pbmc$origin=='PBMC_pre', ]$patient.t <- 'p2.pre'
all.families.pbmc[all.families.pbmc$patient=='p2' & all.families.pbmc$origin=='PBMC_post', ]$patient.t <- 'p2.post'
all.families.pbmc[all.families.pbmc$patient=='p2' & all.families.pbmc$origin=='PBMC_relapse', ]$patient.t <- 'p2.relapse'
all.families.pbmc[all.families.pbmc$patient=='p15' & all.families.pbmc$origin=='PBMC_pre', ]$patient.t <- 'p15.pre'
all.families.pbmc[all.families.pbmc$patient=='p15' & all.families.pbmc$origin=='PBMC_post', ]$patient.t <- 'p15.post'
rownames(all.families.pbmc) <- paste0(all.families.pbmc$patient, '_', rownames(all.families.pbmc))
all.families.pbmc$cell.barcode <- paste0(all.families.pbmc$patient, '_', all.families.pbmc$cell.barcode)

all.families.pbmc <- AddCloneCategory(all.families.pbmc)
dim(all.families.pbmc)

write.table(all.families.pbmc, 'pbmc_combined/tcr_data_raw_pbmc.txt', sep='\t', row.names=F)
all.families.pbmc = read.delim('pbmc_combined/tcr_data_raw_pbmc.txt', sep='\t')
#########################################################################################################
#################################### COMBINE THE TWO DOCUMENTS###########################################
#########################################################################################################

intersecting.rownames = intersect(all.families$cell.barcode, rownames(tils.allmeta)); length(intersecting.rownames)
intersecting.rownames = intersect(all.families.pbmc$cell.barcode, rownames(pbmc.allmeta)); length(intersecting.rownames)

matches.in.tils = all.families %>% subset(cell.barcode %in% rownames(tils.allmeta)); nrow(matches.in.tils)

rownames(pbmc.allmeta) <- pbmc.allmeta$cell.barcode
all.families.in.pbmc = all.families.pbmc %>% subset(cell.barcode %in% rownames(pbmc.allmeta)); nrow(all.families.in.pbmc)

# mit.unique = unique(matches.in.tils$patient.clone); length(mit.unique) 

# tils.allmeta$cell.barcode <- rownames(tils.allmeta)

# match.clons = merge(tils.allmeta[,c('cell.barcode', 'patient.clone.meta', 'TCR.Clone')], 
#                     matches.in.tils[,c('cell.barcode', 'patient.clone', 'TCR.clonotype.ID.revised')],
#                     by='cell.barcode', all.x=TRUE); dim(match.clons)
# View(match.clons)

# matches.in.tils$tcell.clone<-NULL
setdiff(c('cell.barcode', 'patient', 'sample', 'cd4.adt', 'cd8.adt', 'dp.adt', 'dn.adt', 'patient.clone.meta', cite.cols, 'clusters'), colnames(tils.allmeta))
setdiff(c('cell.barcode', 'patient', 'sample', 'cd4.adt', 'cd8.adt', 'dp.adt', 'dn.adt', cite.cols, 'clusters'), colnames(pbmc.allmeta))
table.for.analysis = merge(matches.in.tils, tils.allmeta[,c('cell.barcode', 'patient', 'sample', 'cd4.adt', 'cd8.adt', 'dp.adt', 'dn.adt', 'patient.clone.meta', cite.cols, 'clusters')], all.y=TRUE, by='cell.barcode')
table.for.analysis$patient.x<-NULL; table.for.analysis$patient<-table.for.analysis$patient.y; table.for.analysis$patient.y<-NULL

pbmc.meta.cols.keep = c('cell.barcode', 'patient', 'patient.t', 'sample', 'cd4.adt', 'cd8.adt', 'dp.adt', 'dn.adt', cite.cols, 'clusters')
# matches.cols.keep = setdiff(colnames(all.families.in.pbmc), c('patient', 'patient.t'))
dim(pbmc.allmeta)
table.for.analysis = merge(all.families.in.pbmc, pbmc.allmeta[,pbmc.meta.cols.keep], all.y=TRUE, by=c('cell.barcode', 'patient', 'patient.t'))
dim(table.for.analysis)

table.for.analysis$pbmc.clonotype.patient = paste(table.for.analysis$patient, '-',table.for.analysis$clonotype)
table.for.analysis$pbmc.clonotype.psub = paste(table.for.analysis$patient.t, '-',table.for.analysis$clonotype)

table.for.analysis$til.clonotype.patient = paste(table.for.analysis$patient, '-',table.for.analysis$til.clonotype.family)
table.for.analysis$til.clonotype.psub = paste(table.for.analysis$patient.t, '-',table.for.analysis$til.clonotype.family)


#########################################################################################################
########################################## ADD COUNTS TILS ##############################################
#########################################################################################################
final.table = table.for.analysis

c2 = final.table %>% group_by(patient, patient.clone) %>% tally(name='patient.clonotype.counts'); dim(c2)
nas = c2$patient.clone[is.na(c2$patient.clone)]; length(nas)
c2 = c2[!is.na(c2$patient.clone), ]; dim(c2)
table.for.analysis = merge(table.for.analysis, c2, all.x=TRUE); dim(table.for.analysis)

# 
# c1 = final.table %>% group_by(patient.t, til.clonotype.psub) %>% tally(name='psub.clonotype.counts'); dim(c1)
# nas = grep('- NA', c1$til.clonotype.psub); length(nas)
# c1 = c1[-nas, ]; dim(c1)
# c1 = c1[c1$til.clonotype.psub != 'NA - NA', ]; dim(c1)

# table.for.analysis = merge(table.for.analysis, c1, all.x=TRUE); dim(table.for.analysis)

final.table.2 = table.for.analysis
for (tp in unique(c2$patient)) {
  show(tp)
  sub = c2 %>% subset(patient == tp); show(dim(sub))
  cat('\n')
  colnames(sub) <- c('patient', 'patient.clone', paste0(tp, '.counts'))
  show(head(sub))
  final.table.2 = merge(final.table.2, sub, all.x = TRUE)
}


View(final.table.2[,setdiff(colnames(final.table.2), cite.cols)])
table.for.analysis = final.table.2
# table.for.analysis[1:50,c('patient', 'pbmc.clonotype.psub', 'pbmc.clonotype.patient', 'patient.clonotype.counts', 'psub.clonotype.counts', 'p15.pre.counts', 'p15.post.counts', 'p11.counts', 'p6.counts', 'p2.pre.counts', 'p2.post.counts', 'p2.relapse.counts')]


#########################################################################################################
######################################### ADD COUNTS PBMCS ##############################################
#########################################################################################################
final.table = table.for.analysis

c2 = final.table %>% group_by(patient, pbmc.clonotype.patient) %>% tally(name='patient.clonotype.counts'); dim(c2)
nas = grep('- NA', c2$pbmc.clonotype.patient); length(nas)
c2 = c2[-nas,]; dim(c2)

c1 = final.table %>% group_by(patient.t, pbmc.clonotype.psub) %>% tally(name='psub.clonotype.counts'); dim(c1)
nas = grep('- NA', c1$pbmc.clonotype.psub); length(nas)
c1 = c1[-nas, ]; dim(c1)
c1 = c1[c1$pbmc.clonotype.psub != 'NA - NA', ]; dim(c1)

table.for.analysis = merge(table.for.analysis, c1, all.x=TRUE); dim(table.for.analysis)
table.for.analysis = merge(table.for.analysis, c2, all.x=TRUE); dim(table.for.analysis)

final.table.2 = table.for.analysis
for (tp in unique(c1$patient.t)) {
  show(tp)
  sub = c1 %>% subset(patient.t == tp); show(dim(sub))
  cat('\n')
  colnames(sub) <- c('patient.t', 'pbmc.clonotype.psub', paste0(tp, '.counts'))
  show(head(sub))
  final.table.2 = merge(final.table.2, sub, all.x = TRUE)
}


View(final.table.2[,setdiff(colnames(final.table.2), cite.cols)])
table.for.analysis = final.table.2
table.for.analysis[1:50,c('patient.t', 'pbmc.clonotype.psub', 'pbmc.clonotype.patient', 'patient.clonotype.counts', 'psub.clonotype.counts', 'p15.pre.counts', 'p15.post.counts', 'p11.counts', 'p6.counts', 'p2.pre.counts', 'p2.post.counts', 'p2.relapse.counts')]
#########################################################################################################
######################################### APPLY FUNCTIONS ###############################################
#########################################################################################################


table.for.analysis.majcl = AddMajorityClusterFinal(table.for.analysis, 'patient.clone') # for tils
table.for.analysis.majcl = AddMajorityClusterFinal(table.for.analysis, 'pbmc.clonotype.patient') # for pbmcs
unique(table.for.analysis.majcl$clusters); unique(table.for.analysis.majcl$majority.cluster)
View(table.for.analysis.majcl[,c('pbmc.clonotype.patient', 'pbmc.clonotype.psub', 'patient.clonotype.counts', 'psub.clonotype.counts', 'clusters', 'majority.cluster')])
# final.table = AddCloneCategory(table.for.analysis.majcl)

final.table = SplitChains(table.for.analysis.majcl)
final.table = SplitChains(table.for.analysis)
final.table = SplitChains(final.table)
write.table(final.table, 'tils_combined/tils_cd4_clonotype_families_0515.txt', sep='\t', row.names=F)


final.test = final.table
unique.clonotypes.rough = final.test %>% distinct(patient.clone, category, .keep_all = T); dim(unique.clonotypes.rough)
write.table(unique.clonotypes.rough, 'tils_combined/tils_cd4_clonotype_families_unique_0515.txt', sep='\t', row.names=F)
unique(final.table$clusters)
unique(final.table$majority.cluster)

#########################################################################################################
######################################## combine some columns ###########################################
#########################################################################################################

test.1 = final.table %>% subset(patient=='p15' & patient.clonotype.counts > 1); dim(test.1)
clons = unique(test.1$pbmc.clonotype.patient); length(clons)

for (clon in clons) {
  test.sub1 = test.1 %>% subset(pbmc.clonotype.patient==clon); dim(test.sub1)
  if (length(unique(test.sub1$pbmc.clonotype.psub)) >= 2) {
    subs = sort(unique(test.sub1$pbmc.clonotype.psub))
    print(subs)
    subss = test.sub1 %>% distinct(pbmc.clonotype.patient, pbmc.clonotype.psub, psub.clonotype.counts)
    rownames(subss)<-subss$pbmc.clonotype.psub
    test.1[test.1$pbmc.clonotype.psub==subs[1] & test.1$pbmc.clonotype.patient==clon, ]$p15.pre.counts <- subss[subs[2],]$psub.clonotype.counts
    test.1[test.1$pbmc.clonotype.psub==subs[2] & test.1$pbmc.clonotype.patient==clon, ]$p15.post.counts <- subss[subs[1],]$psub.clonotype.counts
  }
}
View(test.1[test.1$pbmc.clonotype.patient=='p15 - 2',c('patient.t', 'pbmc.clonotype.psub', 'pbmc.clonotype.patient', 'patient.clonotype.counts', 'psub.clonotype.counts', 'p15.pre.counts', 'p15.post.counts', 'p11.counts', 'p6.counts', 'p2.pre.counts', 'p2.post.counts', 'p2.relapse.counts')])

###############

test.2 = final.table %>% subset(patient=='p2' & patient.clonotype.counts > 1); dim(test.2)
clons = unique(test.2$pbmc.clonotype.patient); length(clons)

for (clon in clons) {
  test.sub1 = test.2 %>% subset(pbmc.clonotype.patient==clon); dim(test.sub1)
  if (length(unique(test.sub1$pbmc.clonotype.psub)) >= 2) {
    pres = unique(test.sub1$patient.t)

    if ('p2.post' %in% pres) {
      vals = unique(test.sub1$p2.post.counts); vals = vals[!is.na(vals)]
      test.2[test.2$pbmc.clonotype.patient==clon, ]$p2.post.counts <- vals
    }
    if ('p2.pre' %in% pres) {
      vals = unique(test.sub1$p2.pre.counts); vals = vals[!is.na(vals)]
      test.2[test.2$pbmc.clonotype.patient==clon, ]$p2.pre.counts <- vals
    }
    if ('p2.relapse' %in% pres) {
      vals = unique(test.sub1$p2.relapse.counts); vals = vals[!is.na(vals)]
      test.2[test.2$pbmc.clonotype.patient==clon, ]$p2.relapse.counts <- vals
    }
  }
}
s = test.2[test.2$pbmc.clonotype.patient=='p2 - 6',]; dim(s)
View(test.2[test.2$pbmc.clonotype.patient=='p2 - 6',c('patient.t', 'pbmc.clonotype.psub', 'pbmc.clonotype.patient', 'patient.clonotype.counts', 'psub.clonotype.counts', 'p15.pre.counts', 'p15.post.counts', 'p11.counts', 'p6.counts', 'p2.pre.counts', 'p2.post.counts', 'p2.relapse.counts')])

test.3 = rbind(test.1, test.2)
test.4 = final.table %>% subset(!(cell.barcode %in% test.3$cell.barcode)); dim(test.4)
giacomo.columns = rbind(test.3, test.4); dim(giacomo.columns)

write.table(giacomo.columns, 'pbmc_combined/pbmc_cd8_clonotype_families_mod_0504.txt', sep='\t', row.names=F)

unique.clonotypes.rough = giacomo.columns %>% distinct(pbmc.clonotype.psub, category, .keep_all = T); dim(unique.clonotypes.rough)
write.table(unique.clonotypes.rough, 'pbmc_combined/pbmc_cd8_clonotype_families_mod_unique_0504.txt', sep='\t', row.names=F)

#########################################################################################################
#################################### MAJORITY CLUSTER FUNCTION ##########################################
#########################################################################################################


AddMajorityClusterFinal <- function(complete.table.mult, column_name) {
  # show(column_name)
  # show(dim(complete.table.mult))
  singles = complete.table.mult %>% subset(is.na(patient.clonotype.counts) | patient.clonotype.counts <= 1)
  singles$majority.cluster <- singles$clusters
  
  complete.table.mult = complete.table.mult %>% subset(patient.clonotype.counts > 1); dim(complete.table.mult)
  # show(dim(complete.table.mult))
  # show(dim(singles))
  # show(nrow(complete.table.mult)+nrow(singles))
  # 
  # return(FALSE)
  complete.table = complete.table.mult
  unique.clones = unique(complete.table[,column_name]); cat('number of unique clones:', length(unique.clones), '\n')
  unique.clones = unique.clones[!is.na(unique.clones)]; cat('after removing NAs', length(unique.clones), '\n')
  # na.clones = unique.clones[unique.clones %in% c('p2-NA', 'p6-NA', 'p11.pre-NA', 'p11.rel-NA', 'p15-NA')]; cat('number that are NA:', length(na.clones), '\n')
  complete.table$majority.cluster <- 'None'
  counter = 1
  for (clones in unique.clones) {
    main.cluster = 'None'
    
    clones.sub = complete.table[complete.table[,column_name]==clones & !is.na(complete.table[,column_name]), ]
    cluster.dist = table(clones.sub$clusters)
    num.clones = nrow(clones.sub)
    cluster.dist.frac = cluster.dist / num.clones
    cluster.dist.df = data.frame('frac'=cluster.dist.frac)
    cluster.dist.df$frac.Var1 <- as.character(cluster.dist.df$frac.Var1)

    
    most = top_n(cluster.dist.df, 1, wt=frac.Freq)
    if (nrow(most) == 1) {
      main.cluster = as.character(most$frac.Var1)
    }
    else {
      main.cluster = 'None'
    }
    
    if (counter %% 50 == 0 & !is.na(main.cluster)) { show(counter); 
      show(cluster.dist.df) 
      show(most)
      show(main.cluster) 
      cat('\n')}
    
    complete.table[which(complete.table[,column_name]==clones & !is.na(complete.table[,column_name])), ]$majority.cluster <- main.cluster
    counter = counter + 1
  }
  
  complete.table = rbind(complete.table, singles)
  show(nrow(complete.table))
  return(complete.table)
  
}

SplitChains <- function(clonotypes) {
  show('SPLITTING CHAINS')
  # cat('column names in clonotypes doc:', colnames(clonotypes), '\n')
  clonotypes = clonotypes %>% separate(alpha.1, ',', into=c('TRAV_1', 'TRAJ_1', 'CDR3A_1'))
  clonotypes = clonotypes %>% separate(alpha.2, ',', into=c('TRAV_2', 'TRAJ_2', 'CDR3A_2'))
  clonotypes = clonotypes %>% separate(beta.1, ',', into=c('TRBV_1', 'TRBD_1', 'TRBJ_1', 'CDR3B_1'))
  clonotypes = clonotypes %>% separate(beta.2, ',', into=c('TRBV_2', 'TRBD_2', 'TRBJ_2', 'CDR3B_2'))
  return(clonotypes)
}


AddCloneCategory <- function(all.clonotypes) {
  all.clonotypes$category <- '0'
  show('sub1')
  sub1 = all.clonotypes %>% subset(alpha.1 != 0 & alpha.2 != 0 & beta.1 != 0 & beta.2 == 0)
  sub1$category <- '2.2A1B'
  show('sub2')
  # sub2 = all.clonotypes %>% subset(alpha.1 != 0 & alpha.2 == 0 & beta.1 != 0 & beta.2 != 0)
  # sub2$cateogry <- '3.1A2B'
  show('sub3')
  sub3 = all.clonotypes %>% subset(alpha.1 != 0 & alpha.2 == 0 & beta.1 != 0 & beta.2 == 0)
  sub3$category <- '1.1A1B'
  show('sub4')
  sub4 = all.clonotypes %>% subset((alpha.1 != 0 | alpha.2 != 0) & beta.1 == 0 & beta.2 == 0)
  sub4$category <- '4.A0B'
  show('sub5')
  sub5 = all.clonotypes %>% subset(alpha.1 == 0 & alpha.2 == 0 & (beta.1 != 0 | beta.2 != 0))
  sub5$category <- '5.B0A'
  # show('rest')
  # rest = all.clonotypes %>% subset(is.na(alpha.1)); dim(rest)
  # rest$category <- 'None'
  show('merging all')
  alltog = rbind(sub1, sub3, sub4, sub5)
  # alltog2 = rbind(alltog, rest)
  
  show(nrow(all.clonotypes))
  show(nrow(alltog))
  # show(nrow(alltog2))

  return(alltog)
}















