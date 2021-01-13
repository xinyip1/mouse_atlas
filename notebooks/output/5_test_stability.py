#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
from imp import reload
reload(functions)


# In[4]:


expression  = pd.read_csv('../data/interim/mouse_atlas/mouse_atlas_expression.tsv', sep='\t', index_col=0)
genes       = pd.read_csv('../data/interim/mouse_atlas/mouse_atlas_genes.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('../data/interim/mouse_atlas/mouse_atlas_samples.tsv', sep='\t', index_col=0)
annotations = annotations.rename(columns = {'Platform':'Platform_Category', 'Dataset Name': 'Dataset'})


# In[6]:


# genes = functions.calculate_platform_dependence(blood_data, blood_annotations)

H_index_list, retained_genes_list = functions.resample_clustering(expression, annotations, resample_strategy='jackknife', n_resamples=1, n_clusters_list=[3,4,5,6,7,8])


# In[7]:


H_index_list


# In[11]:


# calculate the mean H_index of each number of cluster k 
for i in H_index_list:
    print( len(i), round(np.mean(i),3), round(min(i),2), round(max(i),2) )

