
setwd('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis')

source('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/code_for_submission/clonotype_grouping.R')

# patient 6
p6.sid = c(2:5)
p6.fl = c("PBMC_sel_TCRs.xlsx", "TIL_sel_1_TCRs.xlsx", "TIL_sel_2_TCRs.xlsx", "TIL_sel_3_TCRs.xlsx", "TIL_sel_4_TCRs.xlsx")
p6.sl = c("PBMC_sel_TCRs", "TIL_sel_1", "TIL_sel_2", "TIL_sel_3", "TIL_sel_4")

pref = 'patient6/data/tcr_data/'
p6.fl <- paste0(pref, p6.fl)
p6.sl <- paste0(pref, p6.sl)


data <- LoadData(p6.fl, p6.sl, 'Patient6')


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
