#!/usr/bin/env python
# coding: utf-8

# ### Visualize the integrated mouse atlas using PCA 

# In[1]:


from sklearn.decomposition import PCA
import pandas as pd
import os, json
import atlas, handler
import plotly.graph_objs as go
import plotly.offline as pyo
pyo.init_notebook_mode()


# In[35]:


atl = atlas.Atlas('../data/interim/mouse_atlas/','mouse_atlas')
atl.runPCA()


# In[41]:


def projectData(atl, key="Cell Type"):
    atl.runPCA()
    coords = atl.coords["filtered"]
    samples = atl.samples("filtered")
    traces = []
    for item in (atl.ordering["filtered"][key]):
        if item != '':
            subset = coords.loc[samples[samples[key]==item].index]
            traces.append(go.Scatter(x=subset['x'], y=subset['y'], 
                                    mode="markers", name="%s (%s)" % (item, len(subset)), hoverinfo="text",
                                    text=[",".join(str(i) for i in samples.loc[sampleId, ['Cell Type', 'Cell Lineage', 'Platform', 'Dataset Name']]) for sampleId in subset.index],
                                    marker=dict(size=7, color=atl.colours["filtered"][key][item])))

    layout = go.Layout(title="PCA", width=1000, height=720, plot_bgcolor='rgb(229, 236, 246)')
    fig = go.Figure(data=traces, layout=layout)
    return fig

# Colour by cell type
fig = projectData(atl, 'Cell Type')
fig.show()
fig.write_image('../output/mouse_atlas_coloured_by_celltype.png')


# In[42]:


# Colour by cell lineage
fig = projectData(atl, 'Cell Lineage')
fig.show()
fig.write_image('../output/mouse_atlas_coloured_by_lineage.png')


# In[43]:


# Colour by platform
fig = projectData(atl, 'Platform')
fig.show()
fig.write_image('../output/mouse_atlas_coloured_by_platform.png')

