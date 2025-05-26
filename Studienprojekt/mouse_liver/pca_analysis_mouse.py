import pandas as pd
import anndata as ad
from pprint import pprint as pp
import PCA_func as pc 
import color_map as cc

path_1 = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt' #add your path to the project
path = path_1+'/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)

'''Full log transformed dataset'''
adataX_df = adata.to_df(layer='log_trasformed')
'''PCA of fulldataset by time as sample'''
def pca_fulldataset_time(x: str):
    adataX_df_time = adataX_df.join(adata.obs['time']).reset_index().drop('index', axis=1).set_index('time').sort_index(ascending=True)
    match x:
        case 'variance':
            pc.variance_ratio(data=adataX_df_time, n=6)
        case 'PCA_2d':
            pc.principal_component_2d_timesample(data=adataX_df_time, label='Time in h', c_map=cc.custom_cmap_func())
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=adataX_df_time, label='Time in h',c_map=cc.custom_cmap_func())
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=adataX_df_time,palette='tab20', hue='time')
'''PCA of fulldataset by study as sample'''
def pca_fulldataset_study(x: str):
    adataX_df_study = adataX_df.join(adata.obs['study']).reset_index().drop('index', axis=1).set_index('study')
    match x:
        case 'PCA_2d':
            pc.principal_component_2d_timesample(data=adataX_df_study, label='Study', c_map=cc.custom_cmap_func('c_34'))
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=adataX_df_study, label='Study',c_map=cc.custom_cmap_func('c_34'))
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=adataX_df_study,palette='tab20', hue='study')
'''PCA of fulldataset by Sequencing Type as sample'''
def pca_fulldataset_seq(x: str):
    adataX_df_st = adataX_df.join(adata.obs['Sequencing Type']).reset_index().drop('index', axis=1).set_index('Sequencing Type')
    match x:
        case 'PCA_2d':
            pc.principal_component_2d_timesample(data=adataX_df_st, label='Sequencing Type', c_map='tab20')
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=adataX_df_st, label='Sequencing Type',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=adataX_df_st,palette='tab20', hue='study')

'''Small reduced dataset also log transformed'''
sr_df = adata.obsm['small_reduction']
'''PCA of small reduced dataset with time as samples'''
def pca_small_time(x: str):
    sr_df_time = sr_df.join(adata.obs['time']).reset_index().drop('index', axis=1).set_index('time').sort_index(ascending=True)
    match x:
        case 'variance':
            pc.variance_ratio(data=sr_df_time, n=6)
        case 'PCA_2d':
            pc.principal_component_2d_timesample(data=sr_df_time, label='Time in h', c_map=cc.custom_cmap_func())
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_time, label='Time in h',c_map=cc.custom_cmap_func())
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=sr_df_time,palette='tab20', hue='time')
'''PCA of small reduced dataset with study as samples'''
def pca_small_study(x: str):
    sr_df_study = sr_df.join(adata.obs['study']).reset_index().drop('index', axis=1).set_index('study')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_study, label='Study',c_map=cc.custom_cmap_func('c_34'))
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=sr_df_study,palette='tab20', hue='study')
'''PCA of small reduced dataset with sex as samples'''
def pca_small_sex(x: str):
    sr_df_sex = sr_df.join(adata.obs['Sex']).reset_index().drop('index', axis=1).set_index('Sex')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_sex, label='Sex',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=sr_df_sex,palette='tab20', hue='Sex')
'''PCA of small reduced dataset with light as samples'''
def pca_small_light(x: str):
    sr_df_light = sr_df.join(adata.obs['Light']).reset_index().drop('index', axis=1).set_index('Light')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_light, label='Light',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=sr_df_light,palette='tab20', hue='Light')
'''PCA of small reduced dataset with age as samples'''
def pca_small_age(x: str):
    sr_df_age = sr_df.join(adata.obs['Age (weeks)']).reset_index().drop('index', axis=1).set_index('Age (weeks)')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_age, label='Age (weeks)',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=sr_df_age,palette='tab20', hue='Age (weeks)')
'''PCA of small reduced dataset with Sequencing Type as samples'''
def pca_small_seq(x: str):
    sr_df_st = sr_df.join(adata.obs['Sequencing Type']).reset_index().drop('index', axis=1).set_index('Sequencing Type')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_st, label='Sequencing Type',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=sr_df_st,palette='tab20', hue='Sequencing Type')
'''PCA of small reduced dataset with Inferred Sequencing Type as samples'''
def pca_small_inf_seq(x: str):
    sr_df_ist = sr_df.join(adata.obs['Inferred Sequencing Type']).reset_index().drop('index', axis=1).set_index('Inferred Sequencing Type')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_ist, label='Sequencing Type',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=sr_df_ist,palette='tab20', hue='Sequencing Type')
'''PCA of small reduced dataset with Note as samples'''
def pca_small_note(x: str):
    sr_df_note = sr_df.join(adata.obs['Note']).reset_index().drop('index', axis=1).set_index('Note')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=sr_df_note, label='Note',c_map='tab20')
        case 'PCA_pari_6':
            pc.pairplot_psa_6(data=sr_df_note,palette='tab20', hue='Note')

'''Strong reduced dataset also log transformed'''
strong_df = adata.obsm['strong_reduction']
'''PCA of strong reduced dataset with time as samples'''
def pca_strong_time(x: str):
    strong_df_time = strong_df.join(adata.obs['time']).reset_index().drop('index', axis=1).set_index('time').sort_index(ascending=True)
    match x:
        case 'variance':
            pc.variance_ratio(data=strong_df_time, n=6)
        case 'PCA_2d':
            pc.principal_component_2d_timesample(data=strong_df_time, label='Time in h', c_map=cc.custom_cmap_func())
        case 'PCA_3d': 
            pc.principal_component_3d_timesample(data=strong_df_time, label='Time in h',c_map=cc.custom_cmap_func())
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_time,palette='tab20', hue='time')
'''PCA of strong reduced dataset with study as samples'''
def pca_strong_study(x: str):
    strong_df_study = strong_df.join(adata.obs['study']).reset_index().drop('index', axis=1).set_index('study')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=strong_df_study, label='Study',c_map=cc.custom_cmap_func('c_34'))
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_study,palette='tab20', hue='study')
'''PCA of strong reduced dataset with sex as samples'''
def pca_strong_sex(x: str):
    strong_df_sex = strong_df.join(adata.obs['Sex']).reset_index().drop('index', axis=1).set_index('Sex')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=strong_df_sex, label='Sex',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_sex,palette='tab20', hue='Sex')
'''PCA of strong reduced dataset with light as samples'''
def pca_strong_light(x: str):
    strong_df_light = strong_df.join(adata.obs['Light']).reset_index().drop('index', axis=1).set_index('Light')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=strong_df_light, label='Light',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_light,palette='tab20', hue='Light')
'''PCA of strong reduced dataset with age as samples'''
def pca_strong_age(x: str):
    strong_df_age = strong_df.join(adata.obs['Age (weeks)']).reset_index().drop('index', axis=1).set_index('Age (weeks)')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=strong_df_age, label='Age (weeks)',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_age,palette='tab20', hue='Age (weeks)')
'''PCA of strong reduced dataset with Sequencing Type as samples'''
def pca_strong_seq(x: str):
    strong_df_st = strong_df.join(adata.obs['Sequencing Type']).reset_index().drop('index', axis=1).set_index('Sequencing Type')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=strong_df_st, label='Sequencing Type',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_st,palette='tab20', hue='Sequencing Type')
'''PCA of strong reduced dataset with Inferred Sequencing Type as samples'''
def pca_strong_inf_seq(x: str):
    strong_df_ist = strong_df.join(adata.obs['Inferred Sequencing Type']).reset_index().drop('index', axis=1).set_index('Inferred Sequencing Type')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=strong_df_ist, label='Sequencing Type',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_ist,palette='tab20', hue='Sequencing Type')
'''PCA of strong reduced dataset with Note as samples'''
def pca_strong_note(x: str):
    strong_df_note = strong_df.join(adata.obs['Note']).reset_index().drop('index', axis=1).set_index('Note')
    match x:
        case 'PCA_3d':
            pc.principal_component_3d_timesample(data=strong_df_note, label='Note',c_map='tab20')
        case 'PCA_pair_6':
            pc.pairplot_psa_6(data=strong_df_note,palette='tab20', hue='Note')



'''Write your code here'''
pca_strong_seq(x='PCA_3d')