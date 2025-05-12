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
'''Open DataSets'''
path_1 = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt' #add your path to the project
path = path_1+'/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)
#cs_df = adata.obsm['centered_study'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
'''Dataset for NMF'''
cs_log1_df = adata.obsm['centered_study_log+1'].join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
cs_log1_df_pos = cs_log1_df+5
'''PCA as baseline'''
# pc.principal_component_3d_timesample(data=cs_log1_df_pos, label='Time in h')
'''PCA with randomized SVD'''
# pc.visualization_of_dec_tools_3d(dec=PCA(svd_solver='randomized'), data=cs_df, label='Time in h')
'''TruncatedSVD'''
# pc.visualization_of_dec_tools_3d(dec=TruncatedSVD(n_components=3), data=cs_df, label='Time in h')
'''FactorAnalysis'''
# pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3), data=cs_df, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3, rotation='varimax'), data=cs_df, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=FactorAnalysis(n_components=3, rotation='quartimax'), data=cs_df, label='Time in h', c_map='twilight')
'''FastICA'''
# pc.visualization_of_dec_tools_3d(dec=FastICA(n_components=3), data=cs_df, label='Time in h')
# pc.visualization_of_dec_tools_3d(dec=FastICA(n_components=3), data=cs_df, label='Time in h', c_map='twilight')
'''NMF'''
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3), data=cs_log1_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='itakura-saito'), data=cs_log1_df_pos, label='Time in h', c_map='twilight')
# pc.visualization_of_dec_tools_3d(dec=NMF(n_components=3, solver='mu', beta_loss='kullback-leibler'), data=cs_log1_df_pos, label='Time in h', c_map='twilight')