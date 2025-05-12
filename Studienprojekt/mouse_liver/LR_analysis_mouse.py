import numpy as np
import pandas as pd
import anndata as ad
from pprint import pprint as pp
import PCA_func as pc 
import color_map as cc
from sklearn.linear_model import LinearRegression

#Import anndata object
path_1 = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH' #add your path to the project
path = path_1+'/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)

'''Small reduced dataset also log transformed''' # uncomment following datasets as necessary
sr_df = adata.obsm['small_reduction']
sr_df_time_setup = sr_df.join(adata.obs['time'])
# sr_df_time = sr_df.join(adata.obs['time']).reset_index().drop('index', axis=1).set_index('time').sort_index(ascending=True)
sr_df_study = sr_df.join(adata.obs['study']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['study', 'time'])
# sr_df_sex = sr_df.join(adata.obs['Sex']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Sex', 'time'])
# sr_df_light = sr_df.join(adata.obs['Light']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Light', 'time'])
# sr_df_age = sr_df.join(adata.obs['Age (weeks)']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Age (weeks)', 'time'])
#sr_df_st = sr_df.join(adata.obs['Sequencing Type']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Sequencing Type', 'time'])
# sr_df_ist = sr_df.join(adata.obs['Inferred Sequencing Type']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Inferred Sequencing Type', 'time'])
# sr_df_note = sr_df.join(adata.obs['Note']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Note', 'time'])
'''PCA of through time centered time as samples'''
# pc.principal_component_3d_timesample(data=pc.residual(sr_df_time, ['time']), label='Time in h')
'''PCA of through study centered time as samples'''
x = pc.residual(sr_df_study, ['study']).reset_index().drop(columns='study').set_index('time').sort_index(ascending=True)
pc.principal_component_3d_timesample(data=x, label='Time in h', c_map='twilight')
'''PCA of through sex centered time as samples'''
# x = pc.residual(sr_df_sex, ['Sex']).reset_index().drop(columns='Sex').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through light centered time as samples'''
# x = pc.residual(sr_df_light, ['Light']).reset_index().drop(columns='Light').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through age centered time as samples'''
# x = pc.residual(sr_df_age, ['Age (weeks)']).reset_index().drop(columns='Age (weeks)').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through Sequencing Type centered time as samples'''
# x = pc.residual(sr_df_st, ['Sequencing Type']).reset_index().drop(columns='Sequencing Type').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through Inferred Sequencing Type centered time as samples'''
# x = pc.residual(sr_df_ist, ['Inferred Sequencing Type']).reset_index().drop(columns='Inferred Sequencing Type').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through Note centered time as samples'''
# x = pc.residual(sr_df_note, ['Note']).reset_index().drop(columns='Note').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')

'''PCA of centered features through linear combination of multiple factors with time as samples'''
sr_df = adata.obsm['small_reduction']
sr_df_time_setup = sr_df.join(adata.obs['time'])
x_0=pc.residual_dummy(data=sr_df, data_dummy_1=adata.obs['Sex'], columns='Sex')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
x_1=pc.residual_dummy(data=x_0, data_dummy_1=adata.obs['Light'], columns='Light')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
x_2=pc.residual_dummy(data=x_1, data_dummy_1=adata.obs['Age (weeks)'], columns='Age (weeks)')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
x_3=pc.residual_dummy(data=x_2, data_dummy_1=adata.obs['Sequencing Type'], columns='Sequencing Type').join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
#x_4=pc.residual_dummy(data=x_3, data_dummy_1=adata.obs['Note'], columns='Note')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
#x_5=pc.residual_dummy(data=x_4, data_dummy_1=adata.obs['Inferred Sequencing Type'], columns='Inferred Sequencing Type')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
#x_6 = pc.residual_dummy(data=x_5, data_dummy_1=adata.obs['study'], columns='study').join(adata.obs['time']).reset_index()#.drop(columns='index').set_index('time').sort_index(ascending=True)
pc.principal_component_3d_timesample(data=x_3, label='Time in h', c_map='twilight')

#Control to see if code works correctly
# x_0=pc.residual_dummy(data=sr_df, data_dummy_1=adata.obs['study'], columns='study').join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x_0, label='Time in h')

'''Strong reduced dataset also log transformed''' # uncomment following datasets as necessary
sr_df = adata.obsm['strong_reduction']
sr_df_time_setup = sr_df.join(adata.obs['time'])
# sr_df_time = sr_df.join(adata.obs['time']).reset_index().drop('index', axis=1).set_index('time').sort_index(ascending=True)
sr_df_study = sr_df.join(adata.obs['study']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['study', 'time'])
# sr_df_sex = sr_df.join(adata.obs['Sex']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Sex', 'time'])
# sr_df_light = sr_df.join(adata.obs['Light']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Light', 'time'])
# sr_df_age = sr_df.join(adata.obs['Age (weeks)']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Age (weeks)', 'time'])
# sr_df_st = sr_df.join(adata.obs['Sequencing Type']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Sequencing Type', 'time'])
# sr_df_ist = sr_df.join(adata.obs['Inferred Sequencing Type']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Inferred Sequencing Type', 'time'])
# sr_df_note = sr_df.join(adata.obs['Note']).join(sr_df_time_setup['time']).reset_index().drop('index', axis=1).set_index(['Note', 'time'])
'''PCA of through time centered time as samples'''
# pc.principal_component_3d_timesample(data=pc.residual(sr_df_time, ['time']), label='Time in h')
'''PCA of through study centered time as samples'''
x = pc.residual(sr_df_study, ['study']).reset_index().drop(columns='study').set_index('time').sort_index(ascending=True)
pc.principal_component_3d_timesample(data=x, label='Time in h', c_map='twilight')
'''PCA of through sex centered time as samples'''
# x = pc.residual(sr_df_sex, ['Sex']).reset_index().drop(columns='Sex').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through light centered time as samples'''
# x = pc.residual(sr_df_light, ['Light']).reset_index().drop(columns='Light').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through age centered time as samples'''
# x = pc.residual(sr_df_age, ['Age (weeks)']).reset_index().drop(columns='Age (weeks)').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through Sequencing Type centered time as samples'''
# x = pc.residual(sr_df_st, ['Sequencing Type']).reset_index().drop(columns='Sequencing Type').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through Inferred Sequencing Type centered time as samples'''
# x = pc.residual(sr_df_ist, ['Inferred Sequencing Type']).reset_index().drop(columns='Inferred Sequencing Type').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')
'''PCA of through Note centered time as samples'''
# x = pc.residual(sr_df_note, ['Note']).reset_index().drop(columns='Note').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x, label='Time in h')

'''PCA of centered features through linear combination of multiple factors with time as samples'''
# sr_df = adata.obsm['strong_reduction']
# sr_df_time_setup = sr_df.join(adata.obs['time'])
# x_0=pc.residual_dummy(data=sr_df, data_dummy_1=adata.obs['Sex'], columns='Sex')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# x_1=pc.residual_dummy(data=x_0, data_dummy_1=adata.obs['Light'], columns='Light')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# x_2=pc.residual_dummy(data=x_1, data_dummy_1=adata.obs['Age (weeks)'], columns='Age (weeks)')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# x_3=pc.residual_dummy(data=x_2, data_dummy_1=adata.obs['Sequencing Type'], columns='Sequencing Type')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# x_4=pc.residual_dummy(data=x_3, data_dummy_1=adata.obs['Note'], columns='Note')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# x_5=pc.residual_dummy(data=x_4, data_dummy_1=adata.obs['Inferred Sequencing Type'], columns='Inferred Sequencing Type')#.join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# x_6 = pc.residual_dummy(data=x_5, data_dummy_1=adata.obs['study'], columns='study').join(adata.obs['time']).reset_index().drop(columns='index').set_index('time').sort_index(ascending=True)
# pc.principal_component_3d_timesample(data=x_6, label='Time in h', c_map='twilight')