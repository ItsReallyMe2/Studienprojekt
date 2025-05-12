The project should be approached in the following way:

First take a look at the supplemental folder, 
there you will find all the necessary datasets related to the study.

Then look at anndata_mouse.py. In this file, 
all datasets were moved to the AnnData object and saved as mice.h5ad.

Next, the datasets were preprocessed in a separate file preprocess_mouse.py. 
And additionally saved again as mice.h5ad.

Next the analysis was subdivided into two files first 
without using linear regression to center the data and after using linear regression, 
the files are called preLR_analysis_mouse.py and LR_analysis_mouse.py respectively.

Additional files like PCA_func.py and color_map.py. 
Contain custom functions and color maps that were created by me and used in the project.
