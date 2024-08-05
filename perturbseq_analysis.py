import os
import time
import copy
import numpy as np
import pandas as pd
import scanpy as sc
import pickle as pkl
import seaborn as sns
import matplotlib.pyplot as plt
from gprofiler import GProfiler
from matplotlib.backends.backend_pdf import PdfPages
from pydeseq2.default_inference import DefaultInference
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.cluster import hierarchy
import spm1d
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from tqdm.notebook import tqdm
import networkx as nx
import sys


########## FILTERING DATA ##########


def filter_data(matrix, var_names, obs_names, mt_names, batch_nums):       # matrix is a dataframe whose head and index are genes and cells, respectively,
    """Filtering data; check comments below"""                             # var_names is a list containing the genes (columns) names (HGNC nomenclature),
    start_time = time.time()                                               # obs_names is a list containing the variants assigned to each cell (row),
    if not os.path.exists('./results'): os.mkdir('./results')              # mt_names is a list containing the mitochondrial genes names,
                                                                           # and batch_nums is a list containing the batches that each cell belongs to

    # Checking parameters                                 
    if matrix.shape[0] != len(obs_names): raise ValueError('Number of rows in matrix does not match the length of variants list')
    if matrix.shape[1] != len(var_names): raise ValueError('Number of columns in matrix does not match the length of genes list')
    if matrix.shape[0] != len(obs_names): raise ValueError('Number of rows in matrix does not match the length of batch annotation list')
    if len(set(matrix.dtypes)) != 1 or 'int64' not in list(set(matrix.dtypes)): raise ValueError('Matrix can only contain integer numbers')

    # Creating metadata
    print('Generating metadata...')
    obs_meta = pd.DataFrame(index = obs_names) # obs_meta = pd.DataFrame(index = range(len(matrix.index)))
    obs_meta['obs_names'] = obs_names
    obs_meta['total_counts'] = list(matrix.sum(axis=1))
    obs_meta['n_genes_counts'] = list(matrix.astype(bool).sum(axis=1))
    obs_meta['batch'] = batch_nums
    mt_found = [mt for mt in mt_names if mt in var_names]
    matrix.columns = var_names
    obs_meta['mt_counts'] = list(matrix[mt_found].sum(axis=1)) # obs_meta['mt_counts'] = pd.Index(matrix[mt_found].sum(axis=1))
    obs_meta['pct_mt_counts'] = round(obs_meta['mt_counts'] / obs_meta['total_counts'] * 100, 2)
    var_meta = pd.DataFrame(index=var_names) # var_meta = pd.DataFrame(index = range(len(matrix.columns)))
    var_meta['var_names'] = var_names
    
    # Setting rules to filter
    adata = sc.AnnData(matrix)
    adata.var = var_meta
    adata.obs = obs_meta
    sc.pl.scatter(adata, x='total_counts', y='n_genes_counts')
    try:
        min_count = int(input('Introduce the minimum number of counts that cells must have: '))
        max_count = int(input('Introduce the maximum number of counts that cells must have: '))
        batch_max = int(input('Introduce the median number of counts that batches must have: '))
    except: print('Not a valid input!')

    # Downsampling cells with an unusually high sequencing depth
    print('Downsampling cells with >%i counts...' % max_count)
    sc.pp.downsample_counts(adata, counts_per_cell=max_count)
    matrix = pd.DataFrame(adata.X)

    # Downsampling batches so they have the same median
    print('Downsampling batches...')
    matrix.index = obs_names # Added
    for batch, batch_median in enumerate(obs_meta.groupby(obs_meta['batch'])['total_counts'].median()):
        if batch_median > batch_max:
            batch_factor = batch_max / batch_median
            matrix.loc[obs_meta['batch'] == batch] = (matrix.loc[obs_meta['batch'] == batch] * batch_factor).astype(int)
    obs_meta['total_counts'] = matrix.sum(axis=1)

    # Creating AnnData
    adata = sc.AnnData(matrix)
    adata.var = var_meta
    adata.obs = obs_meta
    
    # Filtering cells with <7000 counts
    print('Filtering cells with <%i counts...' % min_count)
    adata = adata[adata.obs['total_counts'] > min_count]
    #adata.obs.reset_index(drop=True, inplace=True)
    
    # Filtering cells with >20% mitochondrial counts
    print('Filtering cells with >20% mitochondrial counts...')
    adata = adata[adata.obs['pct_mt_counts'] < 20]
    #adata.obs.reset_index(drop=True, inplace=True)
    
    # Filtering AnnData
    print('Filtering cells with <200 genes')
    sc.pp.filter_cells(adata, min_genes=200)

    # Filtering genes present in <3 cells
    print('Filtering genes present in <5% of cells')
    sc.pp.filter_genes(adata, min_cells=adata.obs.shape[0]*0.05)

    # Normalizing data and filtering lowest variable genes
    print('Normalizing data and filtering lowest variable genes...')
    #adata.var.index = list(adata.var['var_names'])
    adata_nor = copy.deepcopy(adata)
    sc.pp.normalize_total(adata_nor, target_sum=1e4)
    sc.pp.log1p(adata_nor)
    sc.pp.highly_variable_genes(adata_nor, min_mean=0.0125, max_mean=4, min_disp=0.5)
    sc.pl.highest_expr_genes(adata_nor, n_top=20, save='.png')
    sc.pp.scale(adata_nor, max_value=10)  # Calculating z-scores
    adata, adata_nor = adata[:,adata_nor.var['highly_variable']], adata_nor[:,adata_nor.var['highly_variable']]
    #adata.var.reset_index(drop=True, inplace=True), adata_nor.var.reset_index(drop=True, inplace=True)

    # Running Principal Component Analysis
    print('Running Principal Component Analysis...')
    sc.tl.pca(adata_nor, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata_nor, log=True, save='.png')

    # Plotting
    sc.pl.scatter(adata, x='total_counts', y='n_genes_counts', color='pct_mt_counts', save = '.png')
    #sc.pl.scatter(adata, x='total_counts', y='pct_mt_counts', save = '1.png')
    #sc.pl.scatter(adata, x='total_counts', y='n_genes_counts', save = '2.png')

    # Saving
    with open('./results/filter_data.pkl', 'wb') as file: pkl.dump(adata, file)
    with open('./results/normalized_data.pkl', 'wb') as file: pkl.dump(adata_nor, file)
    os.system('mv ./figures/highest_expr_genes.png ./results/filter_data_variable_genes.png')
    os.system('mv ./figures/scatter.png ./results/filter_data_counts.png')
    #os.system('mv ./figures/scatter1.png ./results/filter_data_mt_counts.png')
    #os.system('mv ./figures/scatter2.png ./results/filter_data_genes_counts.png')
    os.system('mv ./figures/pca_variance_ratio.png ./results/filter_data_pca.png')
    os.system('rm -r ./figures')

    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return adata, adata_nor  # both are the AnnData objects filtered, but the second is also normalized and data-transformated


########## HIERARCHIES ANALYSIS ##########


def compare_groups(adata, reference):  # adata is the normalized AnnData object obtained from filter_data function,
                                       # reference is a string with the name of the reference group
    """Calculating Hotelling’s T2 statistic, Pearson score, Spearman value, and L1 linkage between each variant and reference group, deriving
    an empirical null distribution of those scores, and then calculating and plotting the threshold for each method as 5% of desired FDR"""
    start_time = time.time()

    def calculate_metrics(adata_df, adata_pca, reference, variant, obs_names):
        """Calculating metrics"""
        result_list = []
        
        # Retrieving info from selected variants
        reference_df = adata_df[adata.obs[obs_names] == reference]
        variant_df = adata_df[adata.obs[obs_names] == variant]
    
        # Converting into bulkified data
        reference_bulk = np.mean(np.array(reference_df),axis=0)
        variant_bulk = np.mean(np.array(variant_df),axis=0)
    
        # Retrieving PCA info for T2 Hotelling method
        reference_pca = adata_pca[adata.obs[obs_names] == reference]
        variant_pca = adata_pca[adata.obs[obs_names] == variant]
    
        # HotellingT2
        T2 = spm1d.stats.hotellings2(reference_pca, variant_pca)
        result_list.append(T2.z)  # T2 statistic 
    
        # bulk.pearson
        true_value = 1 - pearsonr(reference_bulk, variant_bulk)[0]
        result_list.append(true_value)
        
        # bulk.spearman
        true_value = 1 - spearmanr(reference_bulk, variant_bulk)[0]
        result_list.append(true_value)
                
        # bulk.L1
        true_value = np.sum(np.abs(reference_bulk - variant_bulk)) * 1.0 / reference_df.shape[1]
        result_list.append(true_value)
    
        return result_list

    # Setting data
    adata.obsm['X_pca'] *= -1
    result_df = pd.DataFrame(columns = ['reference', 'variant', 'HotellingT2','bulk.pearson','bulk.spearman','bulk.L1'])
    permuted_df = pd.DataFrame(columns = ['reference', 'variant', 'HotellingT2','bulk.pearson','bulk.spearman','bulk.L1'])
    adata.obs['obs_permuted'] = list(adata.obs['obs_names'])
    adata_df = pd.DataFrame(adata.X)
    adata_df.index = adata.obs.index
    adata_pca = pd.DataFrame(adata.obsm['X_pca']).loc[:,:20]
    adata_pca.index = adata.obs.index
    result_index = 0
    permuted_index = 0

    # Calculating metrics
    print(f'Calculating Hotelling’s T2 statistic, Pearson score, Spearman value, and L1 linkage between each variant and {reference}...')
    variants_list = set(adata.obs['obs_names'])
    with tqdm(total=len(variants_list)) as progress_bar:
        for variant in variants_list:
            if variant != reference:
                result_df.loc[result_index] = [reference, variant] + calculate_metrics(adata_df, adata_pca, reference, variant, 'obs_names')
                result_index += 1
                
                # Permuting data
                for permutation in range(10):
                    adata.obs['obs_permuted'] = list(adata.obs['obs_permuted'].sample(frac=1, random_state=permutation, ignore_index=True))
                    permuted_df.loc[permuted_index] = [variant, 'permuted_' + variant] + calculate_metrics(adata_df, adata_pca, reference, variant, 'obs_permuted')
                    permuted_index += 1
            progress_bar.update(1)

    result_df['dataset'] = 'true'
    permuted_df['dataset'] = 'permuted'
    combo_df = pd.concat([result_df, permuted_df])

    # Calculating FDR
    print('Calculating the threshold for each method as 5% of desired FDR and plotting...')
    min_fdr = 1.0 / combo_df.shape[0]
    for method in ['HotellingT2', 'bulk.pearson', 'bulk.spearman', 'bulk.L1']:
        combo_df = combo_df.sort_values(by=method,ascending=False)
        combo_df['FDR_' + method] = 1
        combo_df = combo_df.reset_index(drop=True)
        for i in range(combo_df.shape[0]):
            value = float(combo_df.loc[i, method])
            passed = combo_df[combo_df[method] >= value]
            true_fraction = np.sum(passed['dataset'] == 'true') / result_df.shape[0]
            perm_fraction = np.sum(passed['dataset'] == 'permuted') / permuted_df.shape[0]
            fdr = max(float(perm_fraction / (true_fraction + perm_fraction)), float(min_fdr))       
            combo_df.loc[i,'FDR_' + method] = float(fdr)
        
        # Plotting
        NUM_BINS=20
        max_val = np.max(combo_df[method])
        mybins=[x * max_val / NUM_BINS for x in range(NUM_BINS)]
        true_scores = combo_df.loc[combo_df['dataset'] == 'true', method]
        perm_scores = combo_df.loc[combo_df['dataset'] == 'permuted', method]
        plt.hist(true_scores, color='red', bins=mybins, label='Real data')
        plt.hist(perm_scores, color='black', alpha=0.5, bins=mybins, weights=[len(true_scores)/len(perm_scores)]*len(perm_scores), label='Permuted data')
        plt.grid(False)
        plt.xlabel(method)
        plt.ylabel('Frequency')

        # Defining threshold
        result_df = combo_df.loc[combo_df['dataset'] == 'true']
        fdr_passed = result_df[result_df['FDR_' + method] <= 0.01].sort_values(by=method,ascending=True)
        thresh = float(list(fdr_passed[method])[0])
        
        plt.axvline(x=thresh,color='black',linestyle='dotted')
        plt.title('FDR 0.05 - Threshold %.2f' % thresh)
        plt.legend()
        plt.savefig('./results/fdr_%s.png' % method, dpi=150)
        plt.show()

    scoring_df = combo_df[combo_df['dataset'] == 'true'].drop(columns=['dataset']).reset_index(drop=True)
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return scoring_df  # Dataframe of metric scores calculated for each variant VS reference group. 


def plot_dendogram(adata, reference, scoring_df, h2_thresh, color_thresh):    # adata is the normalized AnnData object obtained from filter_data function,
    """Hierarchical dendogram based on Spearman scores,"""                    # reference is a string with the name of the reference group,
    """hierarchical clustering based on visual inspection changing"""         # scoring_df is the dataframe obtained from compare_groups() function,
    """the threshold parameter, and plotting the heat map"""                  # h2_thresh is the HotellingT2 threshold obtained in compute_fdr() function,
    start_time = time.time()                                                  # color_thresh is a float number to change as desired the number of clusters in the dendogram

    # Hierarchical clustering
    distance_matrix = 1 - scoring_df['bulk.spearman'].abs().values.reshape(-1, 1)
    Z = hierarchy.linkage(distance_matrix, method='average')
    scoring_df['cluster'] = hierarchy.fcluster(Z, t=color_thresh, criterion='distance')
    scoring_df.to_csv('./results/scoring_df.csv', index=None)
    
    # Plotting
    plt.figure(figsize=(10, 6))
    dn = hierarchy.dendrogram(Z, labels=list(scoring_df['variant']), leaf_rotation=90, color_threshold=color_thresh)  # Generating dendrogram
    plt.title('Dendrogram')
    plt.xlabel('Variants')
    plt.ylabel('Distance')
    
    # Highlighting significant variants
    significant_variants = list(scoring_df[scoring_df['HotellingT2'] < h2_thresh]['variant'])
    ax = plt.gca()
    x_labels = ax.get_xmajorticklabels()
    for label in x_labels:
        if label.get_text() in significant_variants:
            label.set_color('red')  # Change this color as needed
    plt.savefig('./results/scoring_dendogram_variants.png', dpi=150)
    plt.show()

    # Creating dataframe with z-scores of genes across variants
    adata_df = pd.DataFrame(adata.X)
    adata_df.index = adata.obs.index
    variants_list = list(adata.obs['obs_names'].unique())
    variants_list.remove(reference)
    variants_df = pd.DataFrame(index=variants_list, columns=adata.var['var_names'])
    for variant in variants_list:
        variants_df.loc[variant] = list(adata_df[adata.obs['obs_names'] == variant].mean())
    variants_df = variants_df.astype(float)  # Make sure all elements are float numbers
    variants_df = variants_df.loc[scoring_df['variant']]  # Reordering variants to match scoring_df order

    # Plotting cluster heat map
    cluster_list = set(scoring_df['cluster'])
    palette = sns.color_palette("muted", len(cluster_list))
    color_dict = dict(zip(cluster_list, palette))
    color_list = pd.Series(list(scoring_df['cluster']), index=scoring_df['variant']).map(color_dict)
    x = 0.25
    heat_map = sns.clustermap(variants_df.transpose(), row_cluster=True, col_cluster=True, col_linkage=Z, col_colors=color_list,
                   figsize=(20,20), cmap='bwr',vmin=-x,vmax=x,
                   cbar_pos=(0.75, 0.9, 0.15, 0.05), cbar_kws={'orientation':'horizontal', "label": "z-score", 'ticks':[-x,0,x]})
    heat_map.ax_heatmap.set_xlabel('Variants')
    heat_map.ax_heatmap.set_ylabel('Genes')
    plt.savefig('./results/scoring_heat_map.png', dpi=150)
    plt.show()
    
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return scoring_df  # scoring_df now contains a 'cluster' column indicating the group in which the variant is placed


########## LEIDEN ANALYSIS ##########


def leiden_clustering(adata, n_pcs, resolution):   # adata is the normalized AnnData object obtained from previous functions
                                                   # n_pcs is the integer number of Principal Components to use for clustering
                                                   # resolution is a float number to change as desired the number of resulting clusters
    """Leiden clustering of cells and genes, and calculate the variant presence (%) in each leiden group"""
    start_time = time.time()

    # Clustering cells
    print('Clustering cells...')
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.leiden(adata, key_added='leiden', resolution=resolution)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['leiden'], s=1, save = '.png')
    os.system('mv ./figures/umap.png ./results/leiden_cells.png')
    os.system('rm -r ./figures')

    # Calculate variant presence (%) in each leiden group
    obs_names = adata.obs['obs_names']
    result_df = pd.DataFrame(index = list(set(obs_names)), columns = range(len(set(adata.obs['leiden']))))
    for variant in result_df.index:
        cluster_nums = [num for flag, num in zip(obs_names == variant, adata.obs['leiden']) if flag]
        for cluster in set(adata.obs['leiden']):
           result_df.loc[variant, int(cluster)] = round((np.count_nonzero(pd.Index(cluster_nums) == cluster) / len(cluster_nums)) * 100, 2)
    print('Variant presence (%) in each leiden group:')
    display(result_df)
    result_df.to_csv('./results/leiden_cells.csv')
    
    # Transpose and clustering genes
    print('Clustering genes...')
    adata = adata.transpose()
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, key_added='leiden')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['leiden'], s=5, save = '.png')
    os.system('mv ./figures/umap.png ./results/leiden_genes.png')
    os.system('rm -r ./figures')
    print('Number of genes in each leiden group:')
    display(adata.obs['leiden'].value_counts())
    adata = adata.transpose()

    with open('./results/normalized_data.pkl', 'wb') as file:
        pkl.dump(adata, file)
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return adata  # adata is an AnnData object obtained from Leiden clustering


def var_umap(adata):  # adata is the AnnData object obtained from previous functions
    """Plot UMAP by variants"""
    start_time = time.time()
    for variant in set(adata.obs['obs_names']):
        sc.pl.scatter(adata, color='obs_names', groups=[variant, 'WT'], size=40, basis='umap', save = '_%s_WT.png' % variant)
    os.system('mv ./figures ./results/var_umap')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")


########## DIFFERENTIAL EXPRESSION ANALYSIS (DEA) ##########


def diff_analysis(adata, reference):  # adata is the non-normalized AnnData object obtained from the filter_data() function
                                      # reference is a string with the name of the reference group
    """Differential expression analysis (DEA) with the library PyDESeq2, a Python implementation of the DESeq2 method in R"""
    start_time = time.time()
    
    # Defining variables
    result_dict = {}
    inference = DefaultInference()
    if not os.path.exists('./results/diff_analysis'): os.mkdir('./results/diff_analysis')
    
    for variant in set(adata.obs['obs_names']):
        if variant != wt_var:
            # Selecting variants
            adata_temp = adata[adata.obs['obs_names'].isin(['WT', variant])]
            adata_temp.obs.reset_index(drop=True, inplace=True)
            adata_temp.obs['obs_names'].replace(to_replace=['WT', variant], value=['A', 'B'], inplace=True)
        
            # Initializing read counts modeling
            dds = DeseqDataSet(counts=pd.DataFrame(adata_temp.X, columns = adata_temp.var['var_names']), metadata=adata_temp.obs, design_factors='obs_names', refit_cooks=True, inference=inference)
            dds.deseq2()                                     # Fitting dispersions and log fold changes
        
            # Statistical analysis
            stat_res = DeseqStats(dds, inference=inference)  # Compute p-values using Wald tests and adjusted p-values for differential expresion
            stat_res.summary()                               # Running the whole statistical analysis, cooks filtering, and multiple testing adjustement
            result_dict[variant] = stat_res
    
    print('Display results of each variant vs WT doing returned_dict[variant].results_df')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_dict  # result_dict is a dict: keys are variants and values are DEA results with regards to WT


def plot_dea(adata_dict):  # result_dict is the dict obtained from the diff_analysis function
    """Plotting the DEA, and saving info only about LFC and p values"""
    start_time = time.time()
    lfc_df = pd.DataFrame()
    padj_df = pd.DataFrame()
    for variant in adata_dict.keys():

        # Plotting
        stat_res = adata_dict[variant]
        stat_res.plot_MA(s=5)
       #plt.gcf().set_size_inches(8,6)
        plt.gcf().set_dpi(100)
        plt.title(variant + '_WT')
        plt.tight_layout()
        plt.savefig('./results/diff_analysis/' + variant + '_WT.png', dpi=150)
        
        # Saving
        lfc_df[variant] = stat_res.results_df['log2FoldChange']
        padj_df[variant] = stat_res.results_df['padj']
    
    lfc_df = lfc_df.fillna(0).transpose()  # Set NaN values in LFC dataframe to 0
    padj_df = padj_df.fillna(1).transpose()  # Set NaN values in p-values dataframe to 1
    lfc_df.to_csv('./results/diff_analysis/diff_lfc.csv')
    padj_df.to_csv('./results/diff_analysis/diff_padj.csv')
    print('For a higher resolution plot of a particular variant, type returned_dict[variant].plot_MA(s=5)')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return lfc_df, padj_df  # both are dataframes containing LFC and corrected p-values info, respectively


########## NETWORK ANALYSIS ##########


def hotnet_analysis(nodes_df, edges_df, lfc_df, padj_df, scoring_df, wt_cluster):  # nodes_df and edges_df are dataframes containing info about the nodes and edges, respectively, of the constructed network in Cytoscape using STRING
                                                                                   # lfc_df and padj_df are the dataframes obtained from plot_dea() function
                                                                                   # scoring_df is the dataframe obtained from plot_dendogram() function
                                                                                   # wt_cluster is the integer number referring to the Wild type cluster
    """Implementation of Hierarchical HotNet algorithm (https://doi.org/10.1093/bioinformatics/bty613) to find significantly altered subnetworks in each cluster"""

    # Setting environment
    start_time = time.time()
    if wt_cluster not in set(scoring_df['cluster']): raise ValueError('wt_cluster must be an integer number referring to an existing cluster')
    path = './results/hotnet_output'
    if not os.path.exists(f'{path}/data'): os.makedirs(f'{path}/data')
    if not os.path.exists(f'{path}/intermediate/network'): os.makedirs(f'{path}/intermediate/network')
    if not os.path.exists(f'{path}/results'): os.makedirs(f'{path}/results')
    scores_type = ['lfc_scores', 'padj_scores']
    for type in scores_type:
        if not os.path.exists(f'{path}/results/{type}'): os.makedirs(f'{path}/results/{type}')
    num_perm = 10
    cluster_dict = {}
    scoring_df = scoring_df.set_index('variant')
    for cluster in set(scoring_df['cluster']):
        if cluster != wt_cluster: cluster_dict[cluster] = scoring_df.index[scoring_df['cluster'] == cluster].tolist()

    # Generating nodes file
    if not os.path.exists(f'{path}/data/network_nodes.tsv'):
        nodes_df.index = pd.RangeIndex(start=1, stop=len(nodes_df) + 1)
        nodes_df['display name'].to_csv(f'{path}/data/network_nodes.tsv', sep="\t", header=False)

    # Generating edges file
    if not os.path.exists(f'{path}/data/network_edges.tsv'):
        nodes_dict = dict(zip(nodes_df['shared name'], nodes_df['display name']))
        edges_df = edges_df['shared name'].str.split(expand=True)[[0, 2]].replace(nodes_dict)
        nodes_dict = dict(zip(nodes_df['display name'], nodes_df.index))
        edges_df.replace(nodes_dict).to_csv(f'{path}/data/network_edges.tsv', sep="\t", index=False, header=False)

    # Generating scores files
    for cluster, variant_list in cluster_dict.items():
        for variant in variant_list:
            if not os.path.exists(f'{path}/data/lfc_scores_{variant}.tsv'): lfc_df.loc[variant].abs().to_csv(f'{path}/data/lfc_scores_{variant}.tsv', sep="\t", header=False)
            if not os.path.exists(f'{path}/data/padj_scores_{variant}.tsv'): (-np.log10(padj_df.loc[variant])).abs().to_csv(f'{path}/data/padj_scores_{variant}.tsv', sep="\t", header=False)
            if not os.path.exists(f'{path}/intermediate/network_scores_{variant}'): os.makedirs(f'{path}/intermediate/network_scores_{variant}')
            for type in scores_type:
                if not os.path.exists(f'{path}/intermediate/network_scores_{variant}/{type}'): os.makedirs(f'{path}/intermediate/network_scores_{variant}/{type}')

    # Constructing similarity matrices
    print('Constructing similarity matrices...')
    if not os.path.exists(f'{path}/intermediate/network/similarity_matrix.h5') or not os.path.exists(f'{path}/intermediate/network/beta.txt'):
        os.system(f'python ./hotnet/construct_similarity_matrix.py -i {path}/data/network_edges.tsv -o {path}/intermediate/network/similarity_matrix.h5 \
        -bof {path}/intermediate/network/beta.txt')
    
    # Permuting networks
    print('Permuting networks...')
    if not os.path.exists(f'{path}/intermediate/network/index_gene_0.tsv'): os.system(f'cp {path}/data/network_nodes.tsv \
    {path}/intermediate/network/index_gene_0.tsv')
    if not os.path.exists(f'{path}/intermediate/network/edge_list_0.tsv'): os.system(f'cp {path}/data/network_edges.tsv \
    {path}/intermediate/network/edge_list_0.tsv')
    for i in range(1, 5):  # Preserve connectivity of the observed graph
        if not os.path.exists(f'{path}/intermediate/network/edge_list_{i}.tsv'):
            os.system(f'python ./hotnet/permute_network.py -i {path}/intermediate/network/edge_list_0.tsv -s {i} -c \
            -o {path}/intermediate/network/edge_list_{i}.tsv')
    for i in range(5, 9):  # Do not preserve connectivity of the observed graph
        if not os.path.exists(f'{path}/intermediate/network/edge_list_{i}.tsv'):
            os.system(f'python ./hotnet/permute_network.py -i {path}/intermediate/network/edge_list_0.tsv -s {i} \
            -o {path}/intermediate/network/edge_list_{i}.tsv')

    # Permuting scores
    print('Permuting scores...')
    for cluster, variant_list in cluster_dict.items():
        for variant in variant_list:
            for type in scores_type:
                if not os.path.exists(f'{path}/intermediate/network_scores_{variant}/{type}/scores_bins.tsv'):
                    os.system(f'cp {path}/data/{type}_{variant}.tsv {path}/intermediate/network_scores_{variant}/{type}/scores_0.tsv')
                    os.system(f'python ./hotnet/find_permutation_bins.py -gsf {path}/intermediate/network_scores_{variant}/{type}/scores_0.tsv -igf \
                    {path}/data/network_nodes.tsv -elf {path}/data/network_edges.tsv -ms  1000 -o \
                    {path}/intermediate/network_scores_{variant}/{type}/scores_bins.tsv')
                for i in range(1, num_perm + 1):
                    if not os.path.exists(f'{path}/intermediate/network_scores_{variant}/{type}/scores_{i}.tsv'):
                        os.system(f'python ./hotnet/permute_scores.py -i  {path}/intermediate/network_scores_{variant}/{type}/scores_0.tsv \
                        -bf {path}/intermediate/network_scores_{variant}/{type}/scores_bins.tsv -s  {i} -o \
                        {path}/intermediate/network_scores_{variant}/{type}/scores_{i}.tsv')

    # Constructing hierarchies...
    print('Constructing hierarchies...')
    for cluster, variant_list in cluster_dict.items():
        for variant in variant_list:
            for type in scores_type:
                for i in range(num_perm + 1):
                    if not os.path.exists(f'{path}/intermediate/network_scores_{variant}/{type}/hierarchy_edge_list_{i}.tsv') or not os.path.exists(f'{path}/intermediate/network_scores_{variant}/{type}/hierarchy_index_gene_{i}.tsv'):
                        os.system(f'python ./hotnet/construct_hierarchy.py -smf  {path}/intermediate/network/similarity_matrix.h5 -igf \
                        {path}/data/network_nodes.tsv -gsf  {path}/intermediate/network_scores_{variant}/{type}/scores_{i}.tsv -helf \
                        {path}/intermediate/network_scores_{variant}/{type}/hierarchy_edge_list_{i}.tsv -higf \
                        {path}/intermediate/network_scores_{variant}/{type}/hierarchy_index_gene_{i}.tsv')

    # Processing hierarchies...
    print('Processing hierarchies...')
    for cluster, variant_list in cluster_dict.items():
        for variant in variant_list:
            for type in scores_type:
                if not os.path.exists(f'{path}/results/{type}/clusters_network_scores_{variant}.tsv') or not os.path.exists(f'{path}/results/{type}/sizes_network_scores_{variant}.pdf'):
                    os.system(f'python ./hotnet/process_hierarchies.py -oelf {path}/intermediate/network_scores_{variant}/{type}/hierarchy_edge_list_0.tsv \
                    -oigf {path}/intermediate/network_scores_{variant}/{type}/hierarchy_index_gene_0.tsv \
                    -pelf $(for i in `seq {num_perm}`; do echo " {path}/intermediate/network_scores_{variant}/{type}/hierarchy_edge_list_"$i".tsv "; done) \
                    -pigf $(for i in `seq {num_perm}`; do echo " {path}/intermediate/network_scores_{variant}/{type}/hierarchy_index_gene_"$i".tsv "; done)\
                    -lsb  1 -cf {path}/results/{type}/clusters_network_scores_{variant}.tsv -pl network {type}_{variant} \
                    -pf {path}/results/{type}/sizes_network_scores_{variant}.pdf')

    # Performing consensus in each variant...
    print('Performing consensus...')
    for cluster, variant_list in cluster_dict.items():
        igf = ' '.join([f'{path}/data/network_nodes.tsv' for i in range(len(variant_list))])
        elf = ' '.join([f'{path}/data/network_edges.tsv' for i in range(len(variant_list))])
        n = ' '.join(['network' for i in range(len(variant_list))])
        for type in scores_type:
            if not os.path.exists(f'{path}/results/{type}/consensus_nodes_cluster_{cluster}.tsv') or not os.path.exists(f'{path}/results/{type}/consensus_edges_cluster_{cluster}.tsv'):
                cf = ' '.join([f'{path}/results/{type}/clusters_network_scores_{variant}.tsv' for variant in variant_list])
                s = ' '.join([f'{type}_{variant}' for variant in variant_list])
                os.system(f'python ./hotnet/perform_consensus.py -cf {cf} -igf {igf} -elf {elf} -n {n} -s {s} -t 2 \
                -cnf {path}/results/{type}/consensus_nodes_cluster_{cluster}.tsv -cef {path}/results/{type}/consensus_edges_cluster_{cluster}.tsv')

    print(f'Finished! Consensus nodes and edges altered in each variant saved in {path}/results/')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")


def clusnet_analysis(nodes_raw, edges_raw, lfc_df, padj_df, scoring_df, wt_cluster, thresh, perb_gene):  # Parameters are the same as hotnet_analysis function, except:
                                                                                                         # thresh is the fraction (float 0-1) of variants of same cluster a gene must be significant in
                                                                                                         # perb_gene is a string with the name of the perturbed gene
    """Retrieving the intersection of padj and lfc consensus subnetworks, annotating subnetworks with path_gprofiler() function, adding the shortest
    interactor paths between disconnected subnetworks, and finding significant genes appearing in a certain number of variants of the same cluster"""
    start_time = time.time()
    path = './results/hotnet_output'
    scores_type = ['lfc_scores', 'padj_scores']
    scoring_df = scoring_df.set_index('variant')
    interactor_perb = False
    for cluster in set(scoring_df['cluster']):
        if cluster != wt_cluster:
            
            # Saving significant lfc values
            print('Cluster %s' % str(cluster))
            lfc_temp = lfc_df[scoring_df['cluster'] == cluster]
            cond = max(len(lfc_temp.index) * thresh, 2)
            sign_counts = (padj_df < 0.05).loc[lfc_temp.index].sum()
            genes_df =  sign_counts >= cond
            genes_list = list(genes_df[genes_df].index)
            cluster_df = lfc_df.loc[lfc_temp.index, genes_list]
            cluster_df.columns.name = 'genes names'
            pd.DataFrame(cluster_df.mean(), columns = ['significant lfc']).to_csv(f'{path}/results/consensus_lfc_cluster_{cluster}.csv')
    
            # Retrieving intersection of lfc and padj consensus subnetworks
            if not os.path.exists(f'{path}/results/{scores_type[0]}/consensus_edges_cluster_{cluster}.tsv') or not os.path.exists(f'{path}/results/{scores_type[1]}/consensus_edges_cluster_{cluster}.tsv'):
                print('Consensus files from Hotnet analysis not found!')
                sys.exit()
            edges_lfc = pd.read_csv(f'{path}/results/{scores_type[0]}/consensus_edges_cluster_{cluster}.tsv', header=None, sep='\t')
            edges_padj = pd.read_csv(f'{path}/results/{scores_type[1]}/consensus_edges_cluster_{cluster}.tsv', header=None, sep='\t')
            edges_temp = pd.merge(edges_lfc, edges_padj, on=[0, 1])
            if not edges_temp.isin([perb_gene]).any().any(): interactor_perb = True
            edges_consensus = edges_temp.values.tolist()
            
            # Finding and annotating connected components (subnetworks)
            G = nx.Graph()
            G.add_edges_from(edges_consensus)
            subnets = list(nx.connected_components(G))
            subnets_dict = {i: list(component) for i, component in enumerate(subnets)}
            gprofiler_df = path_gprofiler(subnets_dict, 0)
            gprofiler_df.to_csv(f'./results/path_gprofiler_cluster_{cluster}.csv', index=None)
            
            # Creating big interaction network
            nodes_dict = dict(zip(nodes_raw['shared name'], nodes_raw['display name']))
            edges_significant = edges_raw['shared name'].str.split(expand=True)[[0, 2]].replace(nodes_dict)
            edges_significant = list(zip(list(edges_significant[0]), list(edges_significant[2])))
            G = nx.Graph()
            G.add_edges_from(edges_significant)
    
            # Finding shortest paths between subnetworks
            source_subnets = list(subnets_dict.keys())[:-1]
            target_subnets = list(subnets_dict.keys())
            interactor_nodes = []
            for source_key in source_subnets:
                target_subnets.remove(source_key)
                for target_key in target_subnets:
                    shortest_paths, shortest_length = [], float('inf')
                    for gene_1 in subnets_dict[source_key]:
                        for gene_2 in subnets_dict[target_key]:
                            if nx.has_path(G, gene_1, gene_2):
                                path = nx.shortest_path(G, source=gene_1, target=gene_2)
                                path_length = len(path)
                                if path_length < shortest_length:
                                    shortest_paths = [path]
                                    shortest_length = path_length
                                elif path_length == shortest_length:
                                    shortest_paths.append(path)
            
                    # Adding new edges between subnetworks
                    for path in shortest_paths:
                        for i in range(len(path) - 1):
                            if not i == 0 and not path[i] in interactor_nodes: interactor_nodes.append(path[i])
                            if path[i] < path[i + 1] and not [path[i], path[i + 1]] in edges_consensus: edges_consensus.append([path[i], path[i + 1]])
                            elif path[i] > path[i + 1] and not [path[i + 1], path[i]] in edges_consensus: edges_consensus.append([path[i + 1], path[i]])

            # Adding perturb gene if not present in subnetworks
            if interactor_perb and not perb_gene in interactor_nodes:
                interactor_nodes.append(perb_gene)
                edges_perb = list(G.edges(perb_gene))
                for path in edges_perb:
                    path = list(path)
                    if path[0] != perb_gene: index = 0
                    else: index = 1
                    if path[index] in edges_temp or path[index] in interactor_nodes:
                        if path[0] < path[1]: edges_consensus.append(path)
                        else: edges_consensus.append(list(reversed(path)))
            
            # Saving consensus edges with added interactors
            path = './results/hotnet_output'
            edges_consensus = pd.DataFrame(edges_consensus, columns = ['genes names 1', 'genes names 2'])
            edges_consensus['interaction'] = 'pp'
            edges_consensus.to_csv(f'{path}/results/consensus_edges_cluster_{cluster}.csv', index=None)
            interactor_nodes = pd.DataFrame(index=interactor_nodes)
            interactor_nodes['interactor'] = 'yes'
            interactor_nodes.index.name = 'genes names'
            interactor_nodes.to_csv(f'{path}/results/consensus_interactors_cluster_{cluster}.csv')
    
            # Plotting cluster heat map
            x = 2
            heat_map = sns.clustermap(cluster_df.transpose(), row_cluster=True, col_cluster=False, figsize=(20,20),cmap='bwr', vmin=-x,vmax=x,
                                      cbar_pos=(0.75, 0.85, 0.15, 0.05), cbar_kws={'orientation':'horizontal', "label": "LFC value", 'ticks':[-x,0,x]})
            heat_map.ax_heatmap.set_xlabel('Variants')
            heat_map.ax_heatmap.set_ylabel('Genes')
            plt.show()
    
            # Calculating coefficient of variation
            mean_series = cluster_df.mean(axis=0)
            std_series = cluster_df.std(axis=0)
            cv_series = ((std_series / mean_series) * 100).abs().round(decimals=2)
            cv_df = pd.DataFrame(cv_series[cv_series > 30], columns=['CV']).sort_values(by='CV', ascending=False)
            cv_df['Significant counts'] = sign_counts
            print('%i genes of %i genes with CV > 30%%' % (len(cv_df.index), len(cv_series)))
            display(cv_df.transpose())

    print(f'CSV files saved in {path}/results/')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")


########## PATHWAY ENRICHMENT ANALYSIS ##########


def path_gprofiler(genes_dict, verbosity=1):  # genes_dict is a dict whose keys are the cluster numbers, and values are list of genes
                                              # verbosity is a number that if different of 0 additional detail are displayed on screen
    """Pathway enrichment analysis (GO and KEGG) with gprofiler library"""
    start_time = time.time()
                              
    # Running enrichment analysis
    result_cl = ['native', 'name', 'p_value', 'term_size', 'query_size', 'intersection_size', 'effective_domain_size', 'precision', 'recall', 'intersections', 'group']
    result_df = pd.DataFrame(columns=result_cl)
    gp = GProfiler(return_dataframe=True)
    for cluster_num, genes_list in genes_dict.items():
        if verbosity != 0: print(f'Running enrichment analysis in group {cluster_num}...')
        result_temp = gp.profile(organism='hsapiens', query=genes_list, sources=['KEGG', 'GO:BP'], no_evidences=False)
        result_temp['group'] = [cluster_num for i in range(len(result_temp))]
        result_df = pd.concat([result_df, result_temp[result_cl]]).reset_index(drop=True)
    
    if verbosity != 0:
        display(result_df)
        print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_df  # Return dataframe with GO and KEGG annotations of each cluster