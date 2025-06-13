import numpy as np
import pandas as pd
import anndata as ad
from pprint import pprint as pp
import PCA_func as pc 
import umap as um
import cyclops as cc

'''Open DataSets'''
path_1 = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt' #add your path to the project
path = path_1+'/mouse_liver_datasets/mice.h5ad'
adata = ad.read_h5ad(path)
cs_df = adata.obsm['centered_study'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)

cyclops_result = cc.cyclops(cs_df.to_numpy())
cyclops_result.train_model()
cyclops_result.plot_true_vs_predicted_phase(true_phase=adata.obs['time'].to_numpy(), data=cs_df.to_numpy())