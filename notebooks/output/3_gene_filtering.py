#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import functions
get_ipython().run_line_magic('load_ext', 'jupyternotify')


# In[ ]:


expression = pd.read_csv('../data/interim/mouse_atlas/mouse_atlas_expression.tsv', sep='\t', index_col=0)
samples    = pd.read_csv('../data/interim/mouse_atlas/mouse_atlas_samples.tsv', sep='\t', index_col=0)
samples = samples.rename(columns = {'Platform':'Platform_Category'})
print(expression.shape, samples.shape)


# In[ ]:


# calculate platform dependence of each gene 
vp = functions.calculate_platform_dependence(expression, samples)
vp = vp.sort_values(by=['VarFraction'])
get_ipython().run_line_magic('notify', '-m "The cell has finished running"')


# In[ ]:


for threshold in thresholds:
    idx = vp.loc[vp.VarFraction < threshold].index.tolist()
    subset = expression.loc[expression.index.isin(idx)]
    
    pca = sklearn.decomposition.PCA(n_components=10)
    coords = pd.DataFrame(pca.fit_transform(subset.transpose()), index=samples.index, 
                          columns = components)
    coords['platform'] = samples.Platform_Category
    
    for i in components:
        groups = {}
        for i_platform in samples['Platform_Category'].unique():
            sel = samples['Platform_Category']==i_platform
            groups[i_platform] = coords.loc[sel, i]
        args = groups.values()
        KWH.loc[threshold, i] = float(scipy.stats.kruskal(*args)[0])

