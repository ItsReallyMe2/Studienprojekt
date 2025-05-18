import numpy as np
import pandas as pd
import anndata as ad
from pprint import pprint as pp
import PCA_func as pc 
import color_map as cc
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import FactorAnalysis
from sklearn.decomposition import FastICA
from sklearn.decomposition import NMF
from sklearn.manifold import TSNE
import umap as um
from sklearn.cluster import DBSCAN

'''Open DataSets'''
path_1 = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt' #add your path to the project
path = path_1+'/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)
#cs_df_r = adata.obsm['small_reduction'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
#cs_df_r_study = adata.obsm['small_reduction'].join(adata.obs['study']).reset_index().drop(columns='index').set_index('study').sort_index(ascending=True)
cs_df_r_time = adata.obsm['small_reduction'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
#cs_df = adata.obsm['centered_study'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
'''Dataset for NMF'''
# cs_df_pos = cs_df + 6
# cs_log1_df = adata.obsm['centered_study_log+1'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# cs_log1_df_pos = cs_log1_df+5

'''Linear Decompositional tools'''
'''PCA as baseline'''
#pc.principal_component_3d_timesample(data=cs_df_r, label='Time in h')
'''PCA with randomized SVD'''
# pc.visualization_of_dec_tools_3d(dec=PCA(svd_solver='randomized'), data=cs_df, label='Time in h')
'''TruncatedSVD'''
# pc.visualization_of_dec_tools_3d(dec=TruncatedSVD(n_components=3), data=cs_df, label='Time in h')
'''FactorAnalysis'''
# pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3), data=cs_df, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3, rotation='varimax'), data=cs_df, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3, rotation='quartimax'), data=cs_df, label='Time in h', c_map='twilight')
'''NMF with log10(x+1)5'''
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3), data=cs_log1_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu'), data=cs_log1_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='itakura-saito'), data=cs_log1_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='kullback-leibler'), data=cs_log1_df_pos, label='Time in h', c_map='twilight')
'''NMF with log10(x+0.1)+6'''
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3), data=cs_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu'), data=cs_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='itakura-saito'), data=cs_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='kullback-leibler'), data=cs_df_pos, label='Time in h', c_map='twilight')


'''Non Linear Decompositional tools'''
'''t-SNE'''
#pc.visualization_of_dec_tools_3d(dec=TSNE(n_components=3, random_state= 32, perplexity=5), data=cs_df, label='Time in h', c_map='twilight')
'''FastICA'''
# pc.visualization_of_dec_tools_3d(dec=FastICA(n_components=3), data=cs_df, label='Time in h')
# pc.visualization_of_dec_tools_3d(dec=FastICA(n_components=3), data=cs_df, label='Time in h', c_map='twilight')
'''UMAP'''
#pc.visualization_of_dec_tools_3d(dec=um.UMAP(n_components=3), data=cs_df_r_study, label='Time in h', c_map=cc.custom_cmap_func(var_name='c_34'))
umap = um.UMAP(n_components=3)
X_umap = umap.fit_transform(cs_df_r_time)
dbscan = DBSCAN(eps=0.38, n_jobs=-1)
clusters = dbscan.fit_predict(X_umap)
pp(np.unique(clusters))
cs_df_r_time['clusters'] = clusters
pp(cs_df_r_time)
#pc.residual(data=cs_df_r_time, columns='clusters')