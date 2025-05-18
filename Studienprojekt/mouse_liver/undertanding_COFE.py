import COFE.analyse as ana
import COFE.plot as pl
import COFE.scpca as sc
import anndata as ad
from pprint import pprint as pp
import matplotlib.pyplot as plt
path_1 = '/Users/maksimsgolubovics/Python_VScode/Studienprojekt' #add your path to the project
path = path_1+'/mouse_liver_noGH/mice.h5ad'
adata = ad.read_h5ad(path)
df = adata.obsm['centered_study']#.join(adata.obs['time']).reset_index().set_index('time').drop('index', axis=1)
features = df.columns.to_list()
X_test = df.sample(n=10, replace=True, random_state=42)
X_train = df.drop(X_test.index, axis=0).to_numpy()
df_time = adata.obs['time']
true_time_rev = df_time.drop(X_test.index, axis=0)
true_time = df_time.drop(true_time_rev.index, axis=0).to_numpy()
X_test = X_test.to_numpy()

X_train_, X_test_, features_, std_ = ana.preprocess_data(X_train=X_train, X_test=X_test, features=features, feature_dim='col')
dic = ana.cross_validate(X_train=X_train_, s_choices=None ,features=features_)
results = ana.predict_time(X_test=X_test_, cv_results=dic)
del results['features']
pp(results)
pl.plot_circular_ordering(results=results)
plt.show()