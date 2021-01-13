#!/usr/bin/env python
# coding: utf-8

# ## Aim: Collect and combine relevant mouse expression data from the Stemformatics data portal 
# 
# Link: https://www.stemformatics.org/workbench/download_multiple_datasets.
# 
# Stemformatics is an established gene expression data portal containing over 420 public gene expression datasets derived from microarray, RNA sequencing and single cell profiling technologies. It includes curated ‘collections’ of data relevant to cell reprogramming, as well as hematopoiesis and leukaemia.
# 
# ### Samples
# 
# Set the serch field to 'species' and use 'Mus musculus' as search key. 
# 
# ### Processing steps 
# 
# - Sample selection 
# - Combine selected datasets based on platforms
# - Combine all selected datasets

# In[1]:


import pandas as pd
import numpy as np
import atlas
import handler
import requests 


# In[2]:


# inspect the samples metadata
samples = pd.read_csv('/Users/monica/Downloads/export_metadata_samples_v7.2.4.tsv', sep='\t', index_col=2)
samples.head()


# Many of the samples are not healthy blood cells (e.g. iPSCs, AML samples etc.). We will need to select for healthy blood cells from the metadata, and download only the selected samples. 

# ### Sample selction 

# In[3]:


def select_samples(samples):
    '''This function takes the Stemformatics samples metadata and returns samples that are annotated to be blood cells.'''
    pos_selected_id = []
    neg_selected_id = []
    patterns_pos = ['lymp', '[Hh]aem', '[Hh]em', 'HSC','[Mm]ono', '[Bb]-*\ *cell', '[Mm]yelo', 'killer',
                    'NK', '[Mm]eg', '[Bb]aso', '[Nn]eut', '[Ee]os', '[Pp]las', '[Ee]ryt', '\W[Tt]-*\ *cell', 'DC', '[Dd]endri', 
                    'phage', 'macr']
    patterns_neg = ['iPS', 'MSC', 'AML', 'reprogram', 'MAPC', 'KO', 'endothelial', 'LPS', 'mutant', 'Dusp', 'LCMV', 'LSK', 'Chaudhury', 'BLSP',
                    'Bruttger']
    for col in samples.columns:
        l = samples[samples[col].notna()]
        for p in patterns_pos:
            pos_selected_id += l[(l[col].astype(str).str.contains(p) == True)].index.tolist()
        for n in patterns_neg: 
            neg_selected_id += l[(l[col].astype(str).str.contains(n) == True)].index.tolist()
    selected = samples.loc[samples.index.isin(set(pos_selected_id))]
    return selected.loc[~selected.index.isin(set(neg_selected_id))]


# In[4]:


selected_samples = select_samples(samples)
print(selected_samples.shape)
selected_samples.head() # 324 samples are selected


# ### Combine datasets based on platforms

# In[5]:


# add platform information of to samples metadata
datasets = pd.read_csv('/Users/monica/Downloads/export_metadata_datasets_v7.2.4.tsv', sep='\t', index_col=0)
selected_ds = datasets[datasets.index.isin(set(selected_samples.ds_id))]
selected_samples = pd.merge(selected_samples, selected_ds[['description', 'platform']], left_on='ds_id', right_index=True)
selected_samples.columns


# In[6]:


# Inspect the distribution of platforms from which the samples data were generated 
selected_samples.groupby(['platform', 'ds_id']).size()


# In[7]:


# Group selected ds_ids based on platforms
# Note that latforms with not enough dataset representation and small sample size, including Illumina Ref-6, GPL81, MoGene2, 
# and MoEx1 are excluded. 

RNAseq_id             = [7224, 7267, 6655, 6767]
Illu_MouseWG6_id      = [6637, 7291]
Affy_Mouse430_id      = [6498, 6658, 6659, 6756, 6988, 6087, 6108, 6300]
Affy_MoGene1_id       = [6264, 6310, 6313, 6455, 6831, 7131]


# In[8]:


def replace_probes(df, probe_mapping):
    '''
    Input: A microarray expression dataframe and and a probe mapping table.
    Output: The expression dataframe with index changed from probe ids to ensembl gene ids according to the supplied mapping table.
    '''
    probe_mapping.columns = ['ensembl']
    return pd.merge(df, probe_mapping, how='inner', left_index=True, right_index=True).set_index('ensembl')

def find_common_genes(dfs):
    '''
    Input: A list of expression dataframes.
    Output: A list of common genes that are appeared in all dataframes.
    '''
    common_genes = dfs[0].index
    for df in dfs:
        common_genes = common_genes.intersection(df.index)       
    return common_genes

def merge_columns(dfs, common_genes):
    '''
    Input: A list of expression datasts and a list of common genes.
    Output: A combined dataframe of all supplied datasets and keep only the common genes. 
    '''
    matrix = []
    for df in dfs:
        df = df.loc[df.index.isin(common_genes)] # Filter by common genes probed in all datasets 
        m = df.groupby(df.index).first()          # Resolve rows with duplicate indexes by grouping together with the mean value 
        matrix.append(m)
    return pd.concat(matrix, axis=1)

def merge_by_platform(id_list, probe_mapping, platform):
    '''
    This function takes a list of dataset ids and a probe mapping table as inputs, and will load the corresponding datasets into 
    variables formatted as dsxxxx with xxxx indicate the ds_id. If the function is supplied with a probe mapping, it will replace 
    the data index from probe ids to gene ids according to the probe mappintg table.
    '''
    
    serverURL = 'api.stemformatics.org'
    headers = {'Content-type': 'application/json'}
    load_failed = []
    dfs = []
    for ds_id in id_list:
        try:
            result = requests.get('https://%s/expression/%s/raw' % (serverURL, ds_id), headers=headers, verify=False).json()
            df = pd.DataFrame(result['data'], index=result['index'], columns=result['columns'])
            if probe_mapping is not None:
                df = replace_probes(df, probe_mapping) 
            dfs.append(df)
        except ValueError:
            load_failed.append(ds_id)
            continue 
    print('Successfully load {} datasets {}'.format(platform, [i for i in id_list if i not in load_failed]))
    return pd.concat(dfs, axis=1)


# ### Combine all selected datasets

# In[9]:


import warnings
warnings.filterwarnings('ignore')

def combine_s4m_ds():   
    '''
    This function will combine all the selected datasets form the stemformatics portal. 
    '''
    # Load the selected s4m datasets from the data portal api and replace the probe ID with gene ID
    genes_WG6        = pd.read_csv('../data/raw/Stemformatics/probes_MouseWG6_v2.tsv', sep='\t', index_col = 0, header = None)
    genes_Mouse430_2 = pd.read_csv('../data/raw/Stemformatics/probes_Mouse430_2.tsv', sep='\t', index_col = 0, header = None)
    genes_MoGene1    = pd.read_csv('../data/raw/Stemformatics/probes_MoGene1.tsv', sep='\t', index_col = 0, header = None)
    
    
    df_RNAseq    = merge_by_platform(RNAseq_id, None, 'RNAseq')
    df_WG6       = merge_by_platform(Illu_MouseWG6_id, genes_WG6, 'Illumina WG6')
    df_Mouse430  = merge_by_platform(Affy_Mouse430_id, genes_Mouse430_2, 'Affymetrix Mouse430_2')
    df_MoGene1   = merge_by_platform(Affy_MoGene1_id, genes_MoGene1, 'Affymetrix MoGene1')
    
    df_Mouse430.columns = df_Mouse430.columns.str.strip('.CEL')
    df_MoGene1.columns = df_MoGene1.columns.str.strip('.CEL')
    
    dfs = [df_RNAseq, df_WG6, df_Mouse430, df_MoGene1]
    
    common_genes = find_common_genes(dfs)
    combined_df = merge_columns(dfs, common_genes)
    output = combined_df.loc[:,combined_df.columns.isin(selected_samples.index)]
    return atlas.rankTransform(output)

expression_s4m = combine_s4m_ds()


# In[10]:


print(expression_s4m.shape)
expression_s4m.head()


# In[11]:


# Since some of the selected datasets failed to be downloaded from the stemformatics api, we will need to select the samples metadata 
# by the samples that are actually downloaded. 
samples_s4m = selected_samples[selected_samples.index.isin(expression_s4m.columns)]
print(samples_s4m.shape)
samples_s4m.head()


# In[12]:


# Save the combined expression data and samples metadata
expression_s4m.to_csv('../data/interim/mouse_integrate/expression_s4m.tsv', sep='\t')
samples_s4m.to_csv('../data/interim/mouse_integrate/samples_s4m.tsv', sep='\t')


# In[13]:


samples_s4m.shape


# In[14]:


samples_s4m[samples_s4m.dataset_name=='Beerman']

