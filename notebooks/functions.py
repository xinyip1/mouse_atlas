import numpy as np
import pandas as pd
import sklearn
import gc
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sklearn.decomposition
from sklearn.cluster import KMeans, AgglomerativeClustering

def transform_to_percentile(dataframe):

    '''
    Apparently this is properly called the spearman rank

    Parameters:
    ----------

    dataframe
        Dataframe containing expression values, index as variables (genes), columns as samples

    Returns:
    -----------

    transformed_dataframe
        Dataframe with expression as rank (percentile) values

    '''

    transformed_dataframe = (dataframe.shape[0] - dataframe.rank(axis=0, ascending=False, na_option='bottom')+1)/dataframe.shape[0]

    return transformed_dataframe


def calculate_platform_dependence(data, annotations):

    '''
    Calculates the fraction of variance due to platform

    Parameters:
    ----------

    data
        Dataframe containing expression values, index as variables (genes), columns as samples

    annotations
        Dataframe containing metadata, index as samples, requires a column 'Platform_Category'

    Returns:
    ----------

    output_df
        Dataframe contains fraction of variance due to platform

    '''

    import warnings
    from statsmodels.tools.sm_exceptions import ConvergenceWarning
    warnings.simplefilter('ignore', ConvergenceWarning)

    output_df = pd.DataFrame(index=data.index, columns=['VarFraction'])

    data = data.copy().transpose().merge(annotations['Platform_Category'], how='left', left_index=True, right_index=True)

    for i_gene in data.columns.values[:-1]:

        md  = smf.mixedlm("%s ~ Platform_Category" %str(i_gene), data=data, groups = data['Platform_Category'])
        mdf = md.fit()

        output_df.loc[i_gene, 'VarFraction'] = mdf.fittedvalues.std()**2/(mdf.fittedvalues.std()**2+mdf.resid.std()**2)

    return output_df


def resample_clustering(data, annotations, resample_strategy, n_resamples=10, n_clusters_list=[3,4]):

    '''
    This cumbersome function performs either jack-knife or bootstrap resampling,

    Parameters
    ---------

    data
        Expression values, index as genes, columns as samples

    annotations
        Metadata, index as samples, columns as properties
        Must have a 'Platform_Category' column and a 'Dataset' column

    resample_strategy
        Either 'bootstrap' or 'jackknife'

    n_resamples
        Number of bootstrap resampling iterations if resample_stragety is 'bootstrap'

    n_clusters_list
        List of cluster parameters, each value is tested for clustering stability independently

    Returns
    ----------

    results
        Numpy list of arrays, each array contains the H index for each cluster

    retained_genes_list
        Numpy list of arrays, each array contains the retained genes for each iteration of the resampling

    '''
    print("Performing resampling upon data of shape:", data.shape)

    # Initial search for platform dependent genes

    base_platform_dependence = calculate_platform_dependence(data, annotations)

    base_genes  = base_platform_dependence.index.values[base_platform_dependence.VarFraction.values<=0.145]
    base_data   = transform_to_percentile(data.loc[base_genes].copy())
    pca         = sklearn.decomposition.PCA(n_components=3, svd_solver='full')
    base_output = pca.fit_transform(base_data.transpose())

    print("Utlising %d genes as baseline expression data\n" %base_genes.shape[0])

    retained_genes_list = [base_genes]

    resampled_clusters_list = [pd.DataFrame(AgglomerativeClustering(n_clusters=i_clusters).fit_predict(base_output), index=base_data.columns, columns=['Base']) \
                                for i_clusters in n_clusters_list]

    if resample_strategy=='bootstrap':
        iterations = np.arange(n_resamples)
    elif resample_strategy=='jackknife':
        iterations = np.arange(annotations['Dataset'].unique().shape[0])

    #### Do Jacknife
    print("Starting resampling\n")

    for i_iter in iterations:

        if resample_strategy=='jackknife':

            print("Omitting dataset %s" %annotations['Dataset'].unique()[i_iter])
            i_annotations = annotations.copy().loc[annotations['Dataset'] != annotations['Dataset'].unique()[i_iter]]
            i_data        = data[i_annotations.index.values].copy()

        if resample_strategy=='bootstrap':

            print("Bootstrap resampling number %d" %i_iter)
            i_annotations = annotations.copy().sample(frac=1.0, replace=True)
            i_data        = data[i_annotations.index.values].copy()

        i_varPart  = calculate_platform_dependence(transform_to_percentile(i_data), i_annotations)
        i_cut_data = transform_to_percentile(i_data.loc[i_varPart.loc[i_varPart['VarFraction']<=0.145].index.values])
        i_output   = pca.fit_transform(i_cut_data.transpose())

        for i in range(len(n_clusters_list)):

            clusterer       = AgglomerativeClustering(n_clusters=n_clusters_list[i])
            i_clustering    = clusterer.fit_predict(i_output)
            i_clustering_df = pd.DataFrame(index=i_data.columns, columns=[i_iter], data = i_clustering.astype(int))

            resampled_clusters_list[i] = resampled_clusters_list[i].merge(i_clustering_df, how='left', left_index=True, right_index=True)

        gc.collect()

        retained_genes_list.append(i_varPart.loc[i_varPart['VarFraction']<=0.145].index.values)

    results = [calc_H_index(resampled_clusters_list[i]) for i in range(len(n_clusters_list))]

    print('Done\n')

    return results, retained_genes_list

def calc_H_index(dataframe):

    '''
    Calculate H index of from a dataframe of resampled iterations

    Parameters
    ---------

    dataframe
        The first column of the dataframe should be the 'true' reference set
        Agnostic as to the class labels used and they don't have to be consistent across resampling iterations

    Returns
    ---------
        Numpy array containing the H index for each cluster in the true reference set

    '''

    h_index_per_cluster = []

    for i_cluster in dataframe.iloc[:,0].unique():

        samples_in_cluster_i = dataframe.loc[dataframe.iloc[:,0].values == i_cluster].index.values
        max_similarity_list  = []

        for i_resample in range(1, dataframe.shape[1]):

            # Should only compare to samples that were present in the resample (samples may be ommitted in a given iteration of resampling)
            samples_to_compare   = np.intersect1d(samples_in_cluster_i, dataframe.iloc[:, i_resample].dropna().index.values)

            temp_similarity = []

            for j_cluster in dataframe.iloc[:,i_resample].dropna().unique():

                samples_in_cluster_j = dataframe.loc[dataframe.iloc[:,i_resample].values == j_cluster].index.values

                temp_similarity.append(np.intersect1d(samples_in_cluster_j, samples_to_compare).shape[0]/np.union1d(samples_in_cluster_j, samples_to_compare).shape[0])

            max_similarity_list.append(max(temp_similarity))

        h_index = max([h for h in np.arange(0,1.0,0.001) if (np.array(max_similarity_list)>=h).astype(int).sum()/(dataframe.shape[1]-1)>=h])

        h_index_per_cluster.append(h_index)

    return np.array(h_index_per_cluster)
