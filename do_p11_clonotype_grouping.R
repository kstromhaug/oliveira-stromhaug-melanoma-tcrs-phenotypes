
setwd('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis')

source('~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/code_for_submission/clonotype_grouping.R')

# patient 11
p11.sid = c(2:7)
p11.fl = c('new_p11_pbmc_TCRs.xlsx', 'new_p11_til_CD45pos_CD3neg_TCRs.xlsx',
           'new_p11_til_pre_CD45pos_CD3pos_1_TCRs.xlsx', 'new_p11_til_pre_CD45pos_CD3pos_2_TCRs.xlsx', 'new_p11_til_pre_CD45pos_CD3pos_3_TCRs.xlsx',
           'new_p11_til_rel_1_TCRs.xlsx', 'new_p11_til_rel_2_TCRs.xlsx')
p11.sl = c('pbmc', 'til_CD45pos_CD3neg',
           'til_pre_CD45pos_CD3pos_1','til_pre_CD45pos_CD3pos_2', 'til_pre_CD45pos_CD3pos_3',
           'til_rel_1', 'til_rel_2')

pref = 'patient11/data/tcr_data/'
p11.fl <- paste0(pref, p11.fl)
p11.sl <- paste0(pref, p11.sl)


data <- LoadData(p11.fl, p11.sl, 'Patient11')


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


