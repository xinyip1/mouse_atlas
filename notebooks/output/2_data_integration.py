#!/usr/bin/env python
# coding: utf-8

# ## Aim: Integrate gene expression data from three online data portals: Stemformatics, ENCODE and Haemosphere. 
# 
# In previous notebooks, we have combined the gene expression data of selected samples in each individual data portals. Here, we will integrate the three pre-combined expression datasets together into a larger dataset of mouse blood. 
# 

# In[1]:


import pandas as pd
import atlas
import handler
import numpy as np


# ### Merge expression tables

# In[2]:


# Load data 
expression_s4m    = pd.read_csv('../data/interim/mouse_integrate/expression_s4m.tsv', sep='\t', index_col=0)
expression_encode = pd.read_csv('../data/interim/mouse_integrate/expression_encode.tsv', sep='\t', index_col=0)
expression_haemosphere = pd.read_csv('../data/interim/mouse_integrate/expression_Haemosphere.tsv', sep='\t', index_col=0)


# In[3]:


print(expression_s4m.shape, expression_encode.shape, expression_haemosphere.shape)


# In[4]:


dfs = [expression_encode, expression_haemosphere, expression_s4m]
common_genes = handler.find_common_genes(dfs)
mouse_atlas_expression = atlas.rankTransform(handler.merge_columns(dfs, common_genes))

print(mouse_atlas_expression.shape)
mouse_atlas_expression.head()


# ### Standardise sample metadata 
# 
# We have encountered several difficulties during the integration of the metadata that associated with the collected datasets. 
# 
# 1. Different datasets may record different pieces of information to describe each sample. e.g. 'tissue' attribute is recorded in some sample metadata but not all.
# 2. Same piece of information can be recorded in different ways. e.g. cell type information is recorded in a single 'celltype' column in the haemosphere metadata, whereas in s4m metadata, there are 4 columns having information related to the cell type of samples.
# 3. Inconsistent format of contents e.g. macrophage is represented as 'BM mac', 'BM macrophage', 'BM-derived macrophage day 0', 'Bone marrow derived macrophage' in a same column. 
# 4. Different type of information might be stored mixedly under the same attributes. e.g. In 'replicate_group_id' column of the s4m metadata, we might find experiment information, cell type and sort markers information about samples.
# 
# To address these issues:
# 1. we will determin a essential list of attributes to describe samples and unify the naming of these attributes. These essential list of attributes are: 
# 
#        Cell Type; Cell Lineage; Description; Dataset Name; Platform 
# 
# 2. Reannotate the metadata so that the content of each attribute is consistent.

# In[5]:


# Load the reannotated sample metadata 
# For metadata of each data collection: stemformatics, encode and haemosphere, we added two manually annotated columns
# 'cell_lineage_anno' and 'celltype_anno' with consistent content format. 
samples_s4m = pd.read_csv('../data/interim/reannotated/samples_s4m_anno.tsv', sep='\t', index_col=0)
samples_encode = pd.read_csv('../data/interim/reannotated/samples_encode_anno.tsv', sep='\t', index_col=0)
samples_haemosphere = pd.read_csv('../data/interim/reannotated/samples_Haemosphere_anno.tsv', sep='\t', index_col=0)


# slice relevant columns from each metadata table and rename in consistent format

# In[6]:


samples_s4m[:3]


# In[7]:


samples_s4m = samples_s4m[['celltype_anno', 'cell_lineage_anno', 'description', 'platform', 'dataset_name']]
samples_s4m.columns = ['Cell Type', 'Cell Lineage', 'Description', 'Platform', 'Dataset Name']


# In[8]:


samples_haemosphere[:3]


# In[9]:


samples_haemosphere = samples_haemosphere[['celltype_anno', 'cell_lineage_anno', 'description', 'platform', 'dataset_name']]
samples_haemosphere.columns = ['Cell Type', 'Cell Lineage', 'Description', 'Platform', 'Dataset Name']


# In[10]:


samples_encode[:3]


# In[11]:


samples_encode = samples_encode[['celltype_anno', 'cell_lineage_anno', 'Description', 'Platform', 'Project']]
samples_encode.columns = ['Cell Type', 'Cell Lineage', 'Description', 'Platform', 'Dataset Name']


# ### Merge the three metadata tables 

# In[12]:


mouse_atlas_samples = pd.concat([samples_encode, samples_haemosphere, samples_s4m])
print(mouse_atlas_samples.shape)
mouse_atlas_samples.head()


# In[13]:


# remove duplicated samples that are included in multiple data collection
mouse_atlas_samples = mouse_atlas_samples[~mouse_atlas_samples.index.duplicated() & mouse_atlas_samples.index.notna()]
mouse_atlas_expression = mouse_atlas_expression.loc[:,~mouse_atlas_expression.columns.duplicated()]
print(mouse_atlas_samples.shape, mouse_atlas_expression.shape)


# In[14]:


# same the integrated expression table and sample table 
mouse_atlas_expression.to_csv('../data/interim/mouse_atlas/mouse_atlas_expression.tsv', sep='\t')
mouse_atlas_samples.to_csv('../data/interim/mouse_atlas/mouse_atlas_samples.tsv', sep='\t')

