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
from sklearn.cluster import HDBSCAN

'''Open DataSets'''
path_1 = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt' #add your path to the project
path = path_1+'/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)
cs_df = adata.obsm['centered_study'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)


'''Linear Decompositional tools'''
'''PCA as baseline'''
def pca_baseline():
    cs_df_r = adata.obsm['small_reduction'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
    pc.principal_component_3d_timesample(data=cs_df_r, label='Time in h')
'''PCA with randomized SVD'''
def pca_rand_svd():
    pc.visualization_of_dec_tools_3d(dec=PCA(svd_solver='randomized'), data=cs_df, label='Time in h')
'''TruncatedSVD'''
def truncatedSVD():
    pc.visualization_of_dec_tools_3d(dec=TruncatedSVD(n_components=3), data=cs_df, label='Time in h')
'''FactorAnalysis'''
def factorAnalysis(x: str):
    match x:
        case 'no_rotation':
            pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3), data=cs_df, label='Time in h', c_map='twilight')
        case 'varimax':
            pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3, rotation='varimax'), data=cs_df, label='Time in h', c_map='twilight')
        case 'quartimax':
            pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3, rotation='quartimax'), data=cs_df, label='Time in h', c_map='twilight')
'''NMF with log10(x+0.1)+6'''
def nmf(x: str):
    cs_df_pos = cs_df + 6
    match x:
        case 'standard':
            pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3), data=cs_df_pos, label='Time in h', c_map='twilight')
        case 'mu':
            pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu'), data=cs_df_pos, label='Time in h', c_map='twilight')
        case 'mu_itakura-saito':
            pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='itakura-saito'), data=cs_df_pos, label='Time in h', c_map='twilight')
        case 'mu_kullback-leibler':
            pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='kullback-leibler'), data=cs_df_pos, label='Time in h', c_map='twilight')

'''Non Linear Decompositional tools'''
'''t-SNE'''
def tsne():
    pc.visualization_of_dec_tools_3d(dec=TSNE(n_components=3, random_state= 32, perplexity=5), data=cs_df, label='Time in h', c_map='twilight')
'''FastICA'''
def fastICA():
    pc.visualization_of_dec_tools_3d(dec=FastICA(n_components=3), data=cs_df, label='Time in h', c_map='twilight')

'''UMAP'''
def umap_vis():
    pc.visualization_of_dec_tools_3d(dec=um.UMAP(n_components=3, min_dist=0.0, n_neighbors=15, metric='correlation' , n_jobs=-1), data=cs_df, label='Study', c_map='twilight')
'''Batch-Effect reduction with umap and dbscan'''
def umap_hdbscan_batch():
    cs_df_r_time = adata.obsm['small_reduction'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
    umap = um.UMAP(n_components=10, min_dist=0.0)
    X_umap = umap.fit_transform(cs_df_r_time)
    hdbscan = HDBSCAN(n_jobs=-1)
    clusters = hdbscan.fit_predict(X_umap)
    pp(np.unique(clusters))
    cs_df_r_time['clusters'] = clusters.copy()
    pp(cs_df_r_time)
    residual = pc.residual(data=cs_df_r_time, columns=['clusters'])
    pc.principal_component_3d_timesample(data=residual, label='Time in h', c_map='twilight')
