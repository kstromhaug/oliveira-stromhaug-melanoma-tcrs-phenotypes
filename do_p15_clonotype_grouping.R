
setwd('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis')

source('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/code_for_submission/clonotype_grouping.R')

# patient 2
## patient 15
p15.sid = c(1:6)
p15.fl = c("PBMC_post_CD45pos_CD3pos-TCRs.xlsx", "PBMC_pre_CD45pos_CD3pos-TCRs.xlsx", 
           "TIL_CD45pos_CD3pos-1-TCRs.xlsx", "TIL_CD45pos_CD3pos-2-TCRs.xlsx", "TIL_CD45pos_CD3pos-3-TCRs.xlsx", "TIL_CD45pos_CD3pos-4-TCRs.xlsx")
p15.sl = c("PBMC_post_CD45pos_CD3pos", "PBMC_pre_CD45pos_CD3pos",
           "TIL_CD45pos_CD3pos-1", "TIL_CD45pos_CD3pos-2", "TIL_CD45pos_CD3pos-3", "TIL_CD45pos_CD3pos-4")


pref = 'patient15/data/tcr_data/'
p15.fl <- paste0(pref, p15.fl)
p15.sl <- paste0(pref, p15.sl)


data <- LoadData(p15.fl, p15.sl, 'Patient15')


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
tcrs.grouped <- GroupClonotypeFamilies(tcr.groups, clones, empty)
