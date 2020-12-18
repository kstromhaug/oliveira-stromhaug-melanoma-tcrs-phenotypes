import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lm

adata = sc.read_loom('tils_cd8.loom')

obs = adata.obs[['percent_mito', 'patient', 'seurat_clusters', 'nCount_RNA', 'nFeature_RNA', 'til_clonotype_family', 'final_clonotype_family', 'orig_ident']]
var = adata.var
var = var.assign(gene=var.index)
var = var[['gene', 'Selected']]
ldata = anndata.AnnData(X=adata.X, var=var, obs=obs)


sc.pp.normalize_total(ldata, target_sum=1e4)
sc.pp.log1p(ldata)
sc.pp.highly_variable_genes(ldata, min_mean=0.0125, max_mean=3, min_disp=0.5)
ldata = ldata[:, ldata.var.highly_variable]
sc.tl.pca(ldata, svd_solver='arpack')
sc.pp.neighbors(ldata, n_neighbors=40, n_pcs=20)


### DO THE DIFFUSION MAP
root_cell = 'p15_CGTGTCTGTCACTTCC-8' ## tils cd8
root_cell_index = list(ldata.obs.index).index(root_cell)
ldata.uns['iroot'] = root_cell_index

sc.tl.diffmap(ldata, n_comps=15, copy=False)

sc.tl.dpt(ldata, n_branchings=1, n_dcs=10)
sc.tl.dpt(ldata, n_branchings=1, n_dcs=3, min_group_size=10, copy=True) ## <- following the Yost paper
## increase branches?
sc.tl.dpt(ldata, n_branchings=3, n_dcs=3, min_group_size=10, copy=True)
sc.tl.dpt(ldata, n_branchings=1, n_dcs=10, allow_kendall_tau_shift=False, copy=True)

sc.pl.diffmap(ldata, color=['dpt_groups', 'dpt_pseudotime', 'seurat_clusters', 'Mon_State', 'Mon_Pseudotime'])

sc.pl.dpt_groups_pseudotime(ldata)



