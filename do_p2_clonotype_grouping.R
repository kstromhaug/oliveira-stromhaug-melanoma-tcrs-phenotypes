
setwd('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis')

source('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/code_for_submission/clonotype_grouping_old.R')

# patient 2
p2.sid = c(1:8)
p2.fl = c("PBMC_post_CD45pos_CD3pos-TCRs.xlsx", "PBMC_pre_CD45pos_CD3pos-TCRs.xlsx", 
          "PBMC_relapse_CD45pos_CD3pos-TCRs.xlsx", "TIL_CD45neg-TCRs.xlsx", "TIL_CD45pos_CD3neg-TCRs.xlsx",
          "TIL_CD45pos_CD3pos-1-TCRs.xlsx", "TIL_CD45pos_CD3pos-2-TCRs.xlsx", "TIL_CD45pos_CD3pos-3-TCRs.xlsx")
p2.sl = c("PBMC_post_CD45pos_CD3pos", "PBMC_pre_CD45pos_CD3pos",
          "PBMC_relapse_CD45pos_CD3pos", "TIL_CD45neg", "TIL_CD45pos_CD3neg",
          "TIL_CD45pos_CD3pos-1", "TIL_CD45pos_CD3pos-2", "TIL_CD45pos_CD3pos-3")


pref = 'patient2/data/tcr_data/'
p2.fl <- paste0(pref, p2.fl)
p2.sl <- paste0(pref, p2.sl)


data <- LoadData(p2.fl, p2.sl, 'Patient2')


data.split <- SplitData(data)

reordered <- ReorderAlphasBetas(data.split)

prep.output <- PrepTcrData(reordered)
tcr.groups <- prep.output[[1]]
clones <- prep.output[[2]]
empty <- prep.output[[3]]
print(dim(tcr.groups))
print(dim(clones))
print(dim(empty))
print(empty)

##

tcrs.grouped.old <- GroupClonotypeFamilies(tcr.groups, clones, empty)
write.table(tcrs.grouped.old, 'patient2/R_output/Tables/p2_clonotype_families_OLD_20201106.txt', sep='\t')
old.fam.sizes <- tcrs.grouped.old %>% group_by(clonotype) %>% tally(n='family.size')

# tcrs.grouped <- GroupClonotypeFamilies(tcr.groups, clones, empty)
write.table(tcrs.grouped, 'patient2/R_output/Tables/p2_clonotype_families_20201106.txt', sep='\t')
fam.sizes <- tcrs.grouped %>% group_by(clonotype) %>% tally(n='family.size')

# 
# old.tcrgroup = read.xlsx('patient2/R_output/Tables/tils.clonotype.cluster.spreadsheet.xlsx'); dim(old.tcrgroup); colnames(old.tcrgroup)
# 
# old.tcrgroup = read.delim('patient2/R_output/Tables/p2.pbmc.cd8.clonotypes.split.txt', sep='\t'); dim(old.tcrgroup); colnames(old.tcrgroup)
# 
# old.tcrgroup = read.delim('patient2/R_output/til_majority_cluster_spreadsheet_test0421.txt', sep='\t'); dim(old.tcrgroup); colnames(old.tcrgroup)
# 
# old.tcrgroup = read.xlsx('patient2/R_output/Tables/p2.TIL.CD8.clonotype.cluster.spreadsheet.xlsx'); dim(old.tcrgroup); colnames(old.tcrgroup)
# 
# old.tcrgroup = read.delim(''); dim(old.tcrgroup); colnames(old.tcrgroup)


dim(tcrs.grouped.old)
dim(tcrs.grouped)

length(unique(tcrs.grouped.old$clonotype)); length(unique(tcrs.grouped$clonotype))
# length(unique(tcrs.grouped.old$clonotype)); length(unique(tcrs.grouped$clonotype))
table(tcrs.grouped.old$category); table(tcrs.grouped$category)

old.fams <- GetIndividualFamilies(tcrs.grouped.old)
fams <- GetIndividualFamilies(tcrs.grouped)

old.fams.2 <- tcrs.grouped.old[order(tcrs.grouped.old$category), ] %>% distinct(clonotype, .keep_all=T); dim(old.fams.2)
fams.2 <- tcrs.grouped[order(tcrs.grouped$category), ] %>% distinct(clonotype, .keep_all=T); dim(fams.2)
table(old.fams.2$category); table(fams.2$category)

View(old.fams.2 %>% subset(category == '3.1A2B'))
View(fams.2 %>% subset(category == '3.1A2B'))

table(fams.2$category); table(tcrs.grouped$category)


table(old.fams$category); table(old.fams.2$category)



GetIndividualFamilies <- function(families) {
  
  # newdf = data.frame(columns=colnames(families))
  clones <- unique(families$clonotype)
  cat('number of clones:', length(clones), '\n')
  newdf = data.frame(matrix(nrow=length(clones), ncol=ncol(families)))
  colnames(newdf) <- colnames(families)
  
  i = 1
  for (clone in clones) {
    
    if (i %% 500 == 0) { print(i)}
    cl <- families %>% subset(clonotype == clone)
    if ('2.2A1B' %in% cl$category) {
      clsub <- cl %>% subset(category == '2.2A1B')
      # newdf <- rbind(newdf, clsub[1,])
      newdf[i,] <- clsub[1,]
    }
    else if ('3.1A2B' %in% cl$category) {
      clsub <- cl %>% subset(category == '3.1A2B')
      # newdf <- rbind(newdf, clsub[1,])
      newdf[i,] <- clsub[1,]
    }
    else if ('1.1A1B' %in% cl$category) {
      clsub <- cl %>% subset(category == '1.1A1B')
      # newdf <- rbind(newdf, clsub[1,])
      newdf[i,] <- clsub[1,]
    }
    else if ('4.A0B' %in% cl$category) {
      clsub <- cl %>% subset(category == '4.A0B')
      # newdf <- rbind(newdf, clsub[1,])
      newdf[i,] <- clsub[1,]
    }
    else if ('5.B0A' %in% cl$category) {
      clsub <- cl %>% subset(category == '5.B0A')
      # newdf <- rbind(newdf, clsub[1,])
      newdf[i,] <- clsub[1,]
    }
    i = i + 1
  }
  
  return(newdf)
  
}











