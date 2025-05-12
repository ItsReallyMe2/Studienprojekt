import numpy as np
import pandas as pd
import anndata as ad
from pprint import pprint as pp
import matplotlib.pyplot as plt

#Import anndata object
path = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)

'''Reducing the feature space based on the mean and std'''
'''Back to DataFrame'''
adataX_df = adata.to_df(layer='log_trasformed')
'''Figure out the mean cutoff value'''
# bins = np.linspace(-1, 5, 100)
# plt.hist(adataX_df.mean(), bins=bins)
# plt.show() #Huge peak between -1 and -0.89 and a light increase around 0.68
adataX_df_mean_1 = adataX_df.T[adataX_df.mean()>-0.50].T #Decided to use the same average as in the previous case.
adataX_df_mean_2 = adataX_df.T[adataX_df.mean()>-0.50].T
'''Figure out the std cutoff value'''
# bins = np.linspace(0, 3, 200)
# plt.hist(adataX_df_mean_2.std(), bins=bins)
# plt.show()
adataX_df_mean_std_1 = adataX_df_mean_1.T[adataX_df_mean_1.std()>0.22].T
adataX_df_mean_std_2 = adataX_df_mean_2.T[adataX_df_mean_2.std()>0.42].T

'''Adding obsm to anndata'''
adata.obsm['small_reduction'] = adataX_df_mean_std_1 # 12352 transcripts
adata.obsm['strong_reduction'] = adataX_df_mean_std_2 # 2141 transcripts
'''Save anndata object'''
# adata.write('/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/mice.h5ad', compression='gzip')