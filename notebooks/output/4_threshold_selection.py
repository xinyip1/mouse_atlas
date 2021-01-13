#!/usr/bin/env python
# coding: utf-8

# In[25]:


import pandas as pd
import atlas, functions
import numpy as np
import sklearn
import scipy
get_ipython().run_line_magic('load_ext', 'jupyternotify')


# In[49]:


expression = pd.read_csv('../data/interim/mouse_atlas/mouse_atlas_expression.tsv', sep='\t', index_col=0)
samples    = pd.read_csv('../data/interim/mouse_atlas/mouse_atlas_samples.tsv', sep='\t', index_col=0)


# In[51]:


samples = samples.rename(columns = {'Platform':'Platform_Category'})


# In[52]:


# calculate platform dependence of each gene 
vp = functions.calculate_platform_dependence(expression, samples)
vp = vp.sort_values(by=['VarFraction'])
get_ipython().run_line_magic('notify', '-m "The cell has finished running"')


# In[28]:


samples.shape


# In[31]:


expression.shape


# In[57]:


vp[vp.VarFraction<0.2]


# In[53]:


# use kruskal-wallis test to investigate whether samples from different platforms have same distributions 
thresholds = np.round(np.arange(0.04, 1.01, 0.01),2)
components = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
KWH = pd.DataFrame(index=thresholds, columns = components)

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

get_ipython().run_line_magic('notify', '-m "The cell has finished running"')
KWH


# In[54]:


import seaborn as sns
import matplotlib.pyplot as plt
KWH.columns = np.arange(10)+1

plt.figure(figsize=(4,20))
fig = sns.heatmap(KWH.rank(pct=True)[::-1].astype(float), xticklabels=True, yticklabels=True)
fig.set(xlabel='Component', ylabel='Threshold')


# In[ ]:




