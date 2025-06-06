import numpy as np
import pandas as pd
import anndata as ad
from pprint import pprint as pp
import matplotlib.pyplot as plt
import PCA_func as pc

path = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)

'''Reducing the feature space based on the mean and std'''
'''Back to DataFrame'''
adataX_df = adata.to_df(layer='log_trasformed')
adataZ_df = adata.to_df(layer='log+1_trasformed')

'''Figure out the mean and std cutoff value'''
def look_at_dist(x):
    match x:
        case 'mean':
            bins = np.linspace(-1, 5, 100)
            plt.hist(adataZ_df.mean(), bins=bins)
            plt.show()
        case 'std':
            bins = np.linspace(0, 3, 200)
            plt.hist(adataZ_df_mean.std(), bins=bins)
            plt.show()

'''Reduce by mean'''
adataX_df_mean_1 = adataX_df.T[adataX_df.mean()>-0.50].T
adataX_df_mean_2 = adataX_df.T[adataX_df.mean()>-0.50].T
adataZ_df_mean = adataZ_df.T[adataZ_df.mean()>0.15].T

'''Reduce by std'''
adataX_df_mean_std_1 = adataX_df_mean_1.T[adataX_df_mean_1.std()>0.22].T
adataX_df_mean_std_2 = adataX_df_mean_2.T[adataX_df_mean_2.std()>0.42].T
adataZ_df_mean_std = adataZ_df_mean.T[adataZ_df_mean.std()>0.12].T

'''Adding obsm to anndata'''
adata.obsm['small_reduction'] = adataX_df_mean_std_1 # 12352 transcripts
adata.obsm['strong_reduction'] = adataX_df_mean_std_2 # 2141 transcripts
adata.obsm['small_log+1_reducion'] = adataZ_df_mean_std # 12220 transcripts

'''Saving obsm of residual values after centering around study'''
x_0=pc.residual_dummy(data=adata.obsm['small_reduction'], data_dummy_1=adata.obs['study'], columns='study')
x_1=pc.residual_dummy(data=adata.obsm['small_log+1_reducion'], data_dummy_1=adata.obs['study'], columns='study')
adata.obsm['centered_study_log+1'] = x_1
'''Save anndata object'''
# adata.write('/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/mice.h5ad', compression='gzip')