import numpy as np
import pandas as pd
import anndata as ad
from pprint import pprint as pp

'''Assign the dataset path to variables'''
num_reads = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/supplemental/num_reads.by_sample.txt.gz'
tpm = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/supplemental/tpm.by_sample.txt.gz'
sample_metadata = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/supplemental/sample_metadata.txt'
study_metadata = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/supplemental/study_metadata.txt'

'''Open datasets as DataFrame'''
num_reads_df = pd.read_csv(num_reads, compression='gzip', sep='\t')
tpm_df = pd.read_csv(tpm, compression='gzip', sep='\t')
sample_md_df = pd.read_csv(sample_metadata, sep='\t')
study_md_df = pd.read_csv(study_metadata, sep='\t')
reduction_to24h = {24.0 : 0.0,
                 33.0 : 9.0,
                 30.0 : 6.0,
                 36.0 : 12.0,
                 27.0 : 3.0,
                 26.0 : 2.0,
                 34.0 : 10.0,
                 38.0 : 14.0,
                 42.0 : 18.0,
                 32.0 : 8.0,
                 28.0 : 4.0,
                 40.0 : 16.0,
                 44.0 : 20.0,
                 46.0 : 22.0,
                 47.0 : 23.0,
                 65.0 : 17.0,
                 29.0 : 5.0,
                 59.0 : 11.0,
                 41.0 : 17.0,
                 53.0 : 5.0,
                 35.0 : 11.0,
                 64.0 : 16.0,
                 58.0 : 10.0,
                 52.0 : 4.0}
sample_md_df['time'] = sample_md_df['time'].replace(reduction_to24h)

'''AnnData'''
#Prepare just raw expression counts
tpm_df_raw = tpm_df.iloc[:, 2:]
tpm_df_raw.columns = range(tpm_df_raw.shape[1])
tpm_df_raw = tpm_df_raw.T.to_numpy() #to have samples as obs, meaning as rows + transfom to array

#Initializing the AnnData object.. the tpm counts, obs and var 
adata = ad.AnnData(tpm_df_raw)
adata.obs_names = tpm_df.iloc[:, 2:].columns
adata.var_names = tpm_df.set_index('Name').index

#Adding sample_metada to the observation
adata.obs['study'] = sample_md_df.set_index('sample')['study']
adata.obs['time'] = sample_md_df.set_index('sample')['time']
adata.obs['outlier'] = sample_md_df.set_index('sample')['outlier']
#Adding some study_metadata to the observations
ss_join = sample_md_df.set_index('study').join(study_md_df.set_index('Name')).reset_index().set_index('sample')
adata.obs['Sex'] = ss_join['Sex']
adata.obs['Light'] = ss_join['Light']
adata.obs['Age (weeks)'] = ss_join['Age (weeks)']
adata.obs['Sequencing Type'] = ss_join['Sequencing Type']
adata.obs['Inferred Sequencing Type'] = ss_join['Inferred Sequencing Type']
adata.obs['Note'] = ss_join['Note']
#Adding transcript symbol to variables
adata.var['Symbol'] = tpm_df.set_index('Name')['Symbol']

#Adding log10 transformed data
adata.layers['log_trasformed'] = np.log10(adata.X + 0.1)
adata.layers['log+1_trasformed'] = np.log10(adata.X + 1)

#Save anndata object as file
# adata.write('/Users/maksimsgolubovics/Python_VScode/Studienprojekt/mouse_liver_noGH/mice.h5ad', compression='gzip')