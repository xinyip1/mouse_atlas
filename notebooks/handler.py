import pandas as pd
import atlas
from sklearn.decomposition import PCA
import plotly.express as px

def find_common_genes(dfs):
    common_genes = dfs[0].index
    for df in dfs:
        common_genes = common_genes.intersection(df.index)       
    return common_genes

def merge_columns(dfs, common_genes):
    matrix = []
    for df in dfs:
        df = df.loc[df.index.isin(common_genes)] # Filter by common genes probed in all datasets 
        m = df.groupby(df.index).mean()          # Resolve rows with duplicate indexes by grouping together with the mean value 
        matrix.append(m)
    return pd.concat(matrix, axis=1)

def select_samples(samples, patterns = ['[Ll]ymp', '[Hh]aem', '[Hh]em', 'HSC','[Mm]ono', '[Bb]-*\ *cell', '[Mm]yelo', '[Gg]ranu', '[Mm]ac', 
                                        'NK', '[Kk]ill', '[Mm]eg', '[Bb]aso', '[Nn]eut', '[Ee]os', '[Pp]las', '[Ee]ryt', '[Tt]-*\ *cell', 'DC', '[Dd]endri']):
    sid = []
    #colnames = ['sample_type', 'generic_sample_type', 'final_cell_type', 'parental_cell_type']
    colnames = ['celltype', 'description', 'cell_lineage', 'tissue']
    for p in patterns:
        for col in samples.columns:
            l = samples[samples[col].notna()]
            sid += l[l[col].str.contains(p) == True].index.tolist()
    sel_samples = samples.loc[samples.index.isin(sid),]
    print(len(sel_samples), 'out of', len(samples), 'samples are selected.')
    return sel_samples 

def match_reporter(samples, patterns):
    '''
    This function takes in the samples metadata and patterns to be matched, then report the values that are reported as match.
    '''
    colnames = ['sample_type', 'generic_sample_type', 'final_cell_type', 'parental_cell_type']
    out = {}
    for p in patterns: 
        match = []
        for col in samples.columns:
            l = samples[samples[col].notna()]
            is_match = l[col].str.contains(p) == True
            match += l[is_match][col].tolist()
        out[p] = set(match)
    return out

def select_ds(samples):
    dsid = set(select_samples(samples).index)
    print(len(dsid), 'datasets:', str(dsid), 'are selected')
    return dsid

# Replace the expression index from probe id to ensembl id
def replace_probes(df, probe_mapping):
    probe_mapping.columns = ['ensembl']
    return pd.merge(df, probe_mapping, how='inner', left_index=True, right_index=True).set_index('ensembl')

# Load the selected Mouse430_2 datasets as well as the Goodell and Gene-Expression-Commons datasets.
def load_ds(directory, id_list, probe_mapping):
    load_failed = []
    for filename in os.listdir(directory):
        for ds_id in id_list:
            if filename.startswith(str(ds_id)) and filename.endswith(".gct"):  
                df = pd.read_csv(directory+filename, sep='\t',skiprows=2, index_col = 0)
                try: 
                    df = df.drop(['Description'], axis=1)
                    if probe_mapping is not None:
                        df = replace_probes(df, probe_mapping)
                except KeyError:
                    load_failed.append(ds_id)
                    continue
                globals()['ds%s' % ds_id] = df
    print('Successfully load datasets %s' % [i for i in id_list if i not in load_failed])

def find_common_genes(dfs):
    common_genes = dfs[0].index
    for df in dfs:
        common_genes = common_genes.intersection(df.index)       
    return common_genes

def merge_columns(dfs, common_genes):
    matrix = []
    for df in dfs:
        df = df.loc[df.index.isin(common_genes)] # Filter by common genes probed in all datasets 
        m = df.groupby(df.index).mean()          # Resolve rows with duplicate indexes by grouping together with the mean value 
        matrix.append(m)
    return pd.concat(matrix, axis=1)

def runPCA(df, samples):
    df = atlas.rankTransform(df)
    pca = PCA(n_components=3)
    coords = pd.DataFrame(pca.fit_transform(df.values.T), index = samples.index, columns=['x','y','z'])
    coords = pd.merge(coords, samples, left_index=True, right_index=True)
    return coords

def plotPCA(coords, color):
    fig = px.scatter_3d(coords, x='x', y='y', z='z', color=color)
    fig.update_traces(marker=dict(size=3))
    #fig.show()
    return fig

# Example: plotPCA(runPCA(expression_WG6, samples_WG6), 'dataset_name')
