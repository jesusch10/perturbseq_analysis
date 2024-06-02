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
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from matplotlib.backends.backend_pdf import PdfPages
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from pydeseq2.default_inference import DefaultInference
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.cluster import hierarchy
import spm1d
from scipy.stats import pearsonr
from scipy.stats import spearmanr


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
    obs_meta = pd.DataFrame(index = range(len(matrix.index)))
    obs_meta['obs_names'] = obs_names
    obs_meta['total_counts'] = list(matrix.sum(axis=1))
    obs_meta['n_genes_counts'] = list(matrix.astype(bool).sum(axis=1))
    obs_meta['batch'] = batch_nums
    mt_found = [mt for mt in mt_names if mt in list(var_names)]
    matrix.columns = var_names
    obs_meta['mt_counts'] = pd.Index(matrix[mt_found].sum(axis=1))
    obs_meta['pct_mt_counts'] = round(obs_meta['mt_counts'] / obs_meta['total_counts'] * 100, 2)
    var_meta = pd.DataFrame(index = range(len(matrix.columns)))
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
    adata.obs.reset_index(drop=True, inplace=True)
    
    # Filtering cells with >20% mitochondrial counts
    print('Filtering cells with >20% mitochondrial counts...')
    adata = adata[adata.obs['pct_mt_counts'] < 20]
    adata.obs.reset_index(drop=True, inplace=True)
    
    # Filtering AnnData
    print('Filtering cells with <200 genes')
    sc.pp.filter_cells(adata, min_genes=200)

    # Filtering genes present in <3 cells
    print('Filtering genes present in <3 cells')
    sc.pp.filter_genes(adata, min_cells=3)

    # Filtering 80% lowest variable genes
    print('Filtering 80% lowest variable genes')
    adata.var.index = list(adata.var['var_names'])

    adata_copy = copy.deepcopy(adata)
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=4, min_disp=0.5)
    
    sc.pl.highest_expr_genes(adata_copy, n_top=20, save='.png')
    #sc.pp.highly_variable_genes(adata, n_top_genes=int(len(adata.var.index) * 0.2), flavor='seurat_v3', inplace=True)
    adata = adata[:,adata_copy.var['highly_variable']]
    adata.var.reset_index(drop=True, inplace=True)

    # Plotting
    sc.pl.scatter(adata, x='total_counts', y='pct_mt_counts', save = '1.png')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_counts', save = '2.png')

    # Saving
    with open('./results/filter_data.pkl', 'wb') as file:
        pkl.dump(adata, file)
    os.system('mv ./figures/highest_expr_genes.png ./results/filter_data_variable_genes.png')
    os.system('mv ./figures/scatter1.png ./results/filter_data_mt_counts.png')
    os.system('mv ./figures/scatter2.png ./results/filter_data_genes_counts.png')
    os.system('rm -r ./figures')

    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return adata  # adata is the AnnData object filtered


def louvain_cells(adata):  # adata is the AnnData object obtained from filter_data function,
    """Louvain clustering of cells, and calculate the variant presence (%) in each louvain group"""
    start_time = time.time()

    # Louvain clustering
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata, key_added='louvain')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['louvain'], s=1, save = '.png')
    os.system('mv ./figures/umap.png ./results/louvain_cells.png')
    os.system('rm -r ./figures')

    # Calculate variant presence (%) in each louvain group
    obs_names = adata.obs['obs_names']
    result_df = pd.DataFrame(index = list(set(obs_names)), columns = range(len(set(adata.obs['louvain']))))
    for variant in result_df.index:
        cluster_nums = [num for flag, num in zip(obs_names == variant, adata.obs['louvain']) if flag]
        for cluster in set(adata.obs['louvain']):
           result_df.loc[variant, int(cluster)] = round((np.count_nonzero(pd.Index(cluster_nums) == cluster) / len(cluster_nums)) * 100, 2)
    display(result_df)
    result_df.to_csv('./results/louvain_cells.csv')

    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return adata  # adata is an AnnData object obtained from louvain clustering


def var_umap(adata):  # adata is the AnnData object obtained from previous functions
    """Plot UMAP by variants"""
    start_time = time.time()
    for variant in set(adata.obs['obs_names']):
        sc.pl.scatter(adata, color='obs_names', groups=[variant, 'WT'], size=40, basis='umap', save = '_%s_WT.png' % variant)
    os.system('mv ./figures ./results/var_umap')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")


def louvain_genes(adata):  # adata is the AnnData object obtained from previous functions
    """Louvain clustering of genes expressions"""
    start_time = time.time()
    
    # Transpose and clustering
    adata = adata.transpose()
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata, key_added='louvain')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['louvain'], s=5, save = '.png')
    os.system('mv ./figures/umap.png ./results/louvain_genes.png')
    os.system('rm -r ./figures')
    display(adata.var['louvain'].value_counts())

    with open('./results/louvain_data.pkl', 'wb') as file:
        pkl.dump(adata.transpose(), file)
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return adata.transpose()  # adata is an AnnData object obtained from louvain clustering


def path_gprofiler(adata, mart):  # adata is the AnnData object obtained from louvain_genes function,
                                  # and mart is a txt obtained from ensembl mapping (retrieve the genes ID with BioMart)
    """Pathway enrichment analysis with gprofiler library"""
    start_time = time.time()

    # Checking parameters                                 
    if len(adata.var['var_names']) != len(mart): raise ValueError('Number of genes in adata does not match the length of genes ID list')
    
    result_cl = ['louvain', 'source', 'native', 'name', 'p_value', 'significant', 'description', 'term_size', 'query_size', 'intersection_size',
                 'effective_domain_size', 'precision', 'recall', 'query', 'parents', 'intersections', 'evidences']
    result_df = pd.DataFrame(columns=result_cl)
    louvain_serie = adata.var['louvain']
    for louvain_num in range(len(set(louvain_serie))):
        cluster_genes = [id for flag, id in zip(louvain_serie == str(louvain_num), mart) if flag]
        gp = GProfiler(return_dataframe=True)
        result_temp = gp.profile(organism='hsapiens', query=cluster_genes, sources=['KEGG'], no_evidences=False)
        result_temp['louvain'] = [louvain_num for i in range(len(result_temp.index))]
        result_df = pd.concat([result_df, result_temp])
    display(result_df)
    result_df.to_csv('./results/path_gprofiler.csv')

    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_df  # Return dataframe with pathways annotations of each louvain cluster


def path_goatools(adata, entrez):  # adata is the AnnData object obtained from louvain_genes function,
                                   # and entrez is a txt obtained from NCBI mapping (https://pubchem.ncbi.nlm.nih.gov/upload/tools/)
    """Pathway enrichment analysis with goatools library"""
    start_time = time.time()

    # Checking parameters                                 
    if len(adata.var['var_names']) != len(entrez): raise ValueError('Number of genes in adata does not match the length of genes ID list')

    # Loading requirements 
    fin_gene2go = download_ncbi_associations()          # Download GO
    obo_fname = download_go_basic_obo()                 # Download OBO
    obodag = GODag("go-basic.obo")                      # Load GO DAG
    objanno = Gene2GoReader(fin_gene2go, taxids=[9606])
    ns2assoc = objanno.get_ns2assc()                    # Load genes-GO mapping

    # Pathway enrichment analysis
    goeaobj = GOEnrichmentStudyNS(entrez, ns2assoc, obodag, propagate_counts=False, alpha=0.05, methods=['fdr_bh'])

    result_df = pd.DataFrame(columns=['GO', 'term', 'class', 'p', 'p_corr', 'n_genes', 'louvain'])
    louvain_serie = adata.obs['louvain']
    for louvain_num in set(louvain_serie):
        cluster_genes = [int(id) for flag, id in zip(louvain_serie == louvain_num, entrez) if flag and id != '']
        goea_results_all = goeaobj.run_study(cluster_genes)
        goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
        for x in goea_results_sig:
            result_df.loc[len(result_df)] = [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh, cluster_genes, x.ratio_in_study[0],]
    display(result_df)
    result_df.to_csv('./results/path_goatools.csv')

    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_df  # Return dataframe with pathways annotations of each louvain cluster


def diff_analysis(adata):  # adata is the AnnData object obtained from previous functions
    """Differential expression analysis (DEA) with the library PyDESeq2, a Python implementation of the DESeq2 method in R"""
    start_time = time.time()
    
    # Defining variables
    result_dict = {}
    inference = DefaultInference()
    if not os.path.exists('./results/diff_analysis'): os.mkdir('./results/diff_analysis')
    
    for variant in set(adata.obs['obs_names']):
        if variant not in ['WT', 'Null']:
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


def compare_groups(adata, reference):  # adata is the AnnData object obtained from filter_data function,
                                       # reference is a string with the name of the reference group
    """Calculating Hotelling’s T2 statistic, Pearson score, Spearman value, and L1 linkage between each variant and reference group, and deriving an empirical null distribution of those scores"""
    start_time = time.time()
    print('Calculating z-scores...')
    sc.pp.scale(adata, max_value=10)  # Calculating z-scores
    print('Calculating first 50 elements in the Principal Component Space...')
    sc.tl.pca(adata,svd_solver='arpack',n_comps=50)  # Calculating first 50 elements in the Principal Component Space
    adata.obsm['X_pca'] *= -1
    result_df = pd.DataFrame(columns = ['reference', 'variant', 'HotellingT2','bulk.pearson','bulk.spearman','bulk.L1'])
    permuted_df = pd.DataFrame(columns = ['reference', 'variant', 'HotellingT2','bulk.pearson','bulk.spearman','bulk.L1'])
    
    # Permuting data
    adata.obs['obs_permuted'] = list(adata.obs['obs_names'])
    for permutation in range(10):
        print('Permutation', permutation)
        adata.obs['obs_permuted'] = list(adata.obs['obs_permuted'].sample(frac=1, random_state=permutation, ignore_index=True))

        # Calculating metrics
        index = 0
        for variant in set(adata.obs['obs_names']):
            if variant != reference:

                # Retrieving info from selected variants
                adata_df = pd.DataFrame(adata.X)
                adata_df.index = adata.obs.index
                reference_df = adata_df[adata.obs['obs_names'] == reference]
                variant_df = adata_df[adata.obs['obs_names'] == variant]
                reference_pmdf = adata_df[adata.obs['obs_permuted'] == reference]
                variant_pmdf = adata_df[adata.obs['obs_permuted'] == variant]
            
                # Converting into bulkified data
                reference_bulk = np.mean(np.array(reference_df),axis=0)
                variant_bulk = np.mean(np.array(variant_df),axis=0)
                reference_pmbulk = np.mean(np.array(reference_pmdf),axis=0)
                variant_pmbulk = np.mean(np.array(variant_pmdf),axis=0)
            
                # Retrieving PCA info for T2 Hotelling method
                adata_pca = pd.DataFrame(adata.obsm['X_pca']).loc[:,:20]
                adata_pca.index = adata.obs.index
                reference_pca = adata_pca[adata.obs['obs_names'] == reference]
                variant_pca = adata_pca[adata.obs['obs_names'] == variant]
                reference_pmpca = adata_pca[adata.obs['obs_permuted'] == reference]
                variant_pmpca = adata_pca[adata.obs['obs_permuted'] == variant]
            
                result_list = [reference, variant]
                permuted_list = [variant, 'permuted_' + variant]
            
                # HotellingT2
                T2 = spm1d.stats.hotellings2(reference_pca, variant_pca)
                result_list.append(T2.z)  # T2 statistic 
                T2 = spm1d.stats.hotellings2(reference_pmpca, variant_pmpca)
                permuted_list.append(T2.z)
            
                # bulk.pearson
                true_value = 1 - pearsonr(reference_bulk, variant_bulk)[0]
                result_list.append(true_value)
                permuted_value = 1 - pearsonr(reference_pmbulk, variant_pmbulk)[0]
                permuted_list.append(permuted_value)
                
                # bulk.spearman
                true_value = 1 - spearmanr(reference_bulk, variant_bulk)[0]
                result_list.append(true_value)
                permuted_value = 1 - spearmanr(reference_pmbulk, variant_pmbulk)[0]
                permuted_list.append(permuted_value)
                        
                # bulk.L1
                true_value = np.sum(np.abs(reference_bulk - variant_bulk)) * 1.0 / reference_df.shape[1]
                result_list.append(true_value)
                permuted_value = np.sum(np.abs(reference_pmbulk - variant_pmbulk)) * 1.0 / reference_df.shape[1]
                permuted_list.append(permuted_value)
                
                result_df.loc[index], permuted_df.loc[index] = result_list, permuted_list
                index += 1

    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return(adata, result_df, permuted_df)  # adata is an AnnData object with z-scores, and the rest are both dataframes of metric scores calculated for each variant VS reference group, and as null distribution


def compute_fdr(result_df, permuted_df):  # both arguments are dataframes obtained from compare_groups() function
    """Calculating the threshold for each method as 5% of desired FDR"""
    start_time = time.time()
    result_df['dataset'] = 'true'
    permuted_df['dataset'] = 'permuted'
    combo_df = pd.concat([result_df, permuted_df])
    min_fdr = 1.0 / combo_df.shape[0]
    for method in ['HotellingT2', 'bulk.pearson', 'bulk.spearman', 'bulk.L1']:
        combo_df = combo_df.sort_values(by=method,ascending=False)
        combo_df['FDR_' + method] = 1
        combo_df = combo_df.reset_index(drop=True)

        # Calculating FDR
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
        plt.title('FDR 0.01 - Threshold %.2f' % thresh)
        plt.legend()
        plt.savefig('./results/fdr_%s.png' % method, dpi=150)
        plt.show()

    result_df = combo_df[combo_df['dataset'] == 'true'].drop(columns=['dataset']).reset_index(drop=True)
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_df


def plot_dendogram(result_df, h2_thresh, color_thresh):  # result_df is the dataframe obtained from either compare_groups() or compute_fdr() functions
                                                         # h2_thresh is the HotellingT2 threshold obtained in compute_fdr() function
                                                         # color_thresh is a float number to change as desired the number of cluster in the dendogram
    """Hierarchical dendogram based on Pearson scores, and hierarchical clustering based on visual inspection changing the threshold parameter"""
    start_time = time.time()
    distance_matrix = 1 - result_df['bulk.pearson'].abs().values.reshape(-1, 1)
    Z = hierarchy.linkage(distance_matrix, method='average')  # Hierarchical clustering
    
    # Plotting
    plt.figure(figsize=(10, 6))
    dn = hierarchy.dendrogram(Z, labels=list(result_df['variant']), leaf_rotation=90, color_threshold=color_thresh)  # Generating dendrogram
    plt.title('Dendrogram')
    plt.xlabel('Variants')
    plt.ylabel('Distance')
    
    # Highlighting significant variants
    significant_variants = list(result_df[result_df['HotellingT2'] < h2_thresh]['variant'])
    ax = plt.gca()
    x_labels = ax.get_xmajorticklabels()
    for label in x_labels:
        if label.get_text() in significant_variants:
            label.set_color('red')  # Change this color as needed
    
    plt.savefig('./results/scoring_variants.png', dpi=150)
    plt.show()
    result_df['cluster'] = hierarchy.fcluster(Z, t=color_thresh, criterion='distance')
    result_df.to_csv('./results/scoring_variants.csv')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_df  # result_df now contains a 'cluster' column indicating the group in which the variant is placed


def lfc_cluster(lfc_df, padj_df, scoring_df, thresh, perb_gene):  # lfc_df is the dataframe with LFC info obtained from plot_dea() function
                                                                  # padj_df is the dataframe with corrected p-values obtained from plot_dea() function
                                                                  # scoring_df is the dataframe with cluster info obtained from plot_dendogram() function
                                                                  # thresh is the fraction (float 0-1) of variants of same cluster a gene must be significant in
                                                                  # per_genes is a string with the name of the perturbed gene
    """Creating CSV for Cytoscape with LFC info from significant genes appearing in all variants of the same cluster"""
    start_time = time.time()
    scoring_df = scoring_df.set_index('variant')
    for cluster in set(scoring_df['cluster']):
    
        # Getting cluster info 
        print('Cluster %s' % str(cluster))
        cluster_df = lfc_df[scoring_df['cluster'] == cluster]
        cond = max(len(cluster_df.index) * thresh, 2)
        genes_df = (padj_df < 0.05).loc[cluster_df.index].sum() >= cond
        genes_list = list(genes_df[genes_df].index)
        if len(genes_list) != 0: genes_list.append(perb_gene)
        final_df = lfc_df.loc[cluster_df.index, genes_list]
    
        # Plotting data
        if final_df.shape[1] == 0: print('No significant genes found across %.2f of variants. Maybe cluster %s contains variants similar to WT?' % (thresh, str(cluster)))
        else:
            condition = '0'
            while condition == '0':
                plt.figure(figsize=(10, 6))
                font_size = min(max(10 - (final_df.shape[0] // 10) - (final_df.shape[1] // 10), 4), 10)
                sns.heatmap(final_df, annot=True, cmap='coolwarm', center=0, annot_kws={"size": font_size})
                plt.title('LFC Heatmap of significantly expressed genes appearing in all variants')
                plt.xlabel('Genes')
                plt.ylabel('Variants')
                plt.xticks(fontsize=font_size * 1.2)
                plt.yticks(fontsize=font_size * 1.2)
                plt.show()
        
                # Filtering variants
                condition = str(input('Introduce 0 if you want to filter some variants:'))
                if condition == '0':
                    variants_list = input('Introduce variants names to filter separated by commas:').split(',')
                    if len(variants_list) != 0:
                        for variant in variants_list: 
                            if variant not in final_df.index:
                                print(variant, 'variant not found')
                                variants_list.remove(variant)
                        final_df = final_df.loc[final_df.index.drop(variants_list)]
    
            final_df.columns.name = 'genes_names'
            pd.DataFrame(final_df.mean(), columns = ['lfc_cluster%s' % str(cluster)]).to_csv('./results/lfc_cluster%s.csv' % str(cluster))
    print("Execution time:", round(time.time() - start_time, 3), "seconds")


def hotnet_analysis(nodes_df, edges_df, lfc_df, scoring_df, wt_cluster):  # nodes_df and edges_df are dataframes containing info about the nodes and edges, respectively, of the constructed network in Cytoscape using STRING
                                                                          # lfc_df and scoring_df are dataframes obtained from plot_dea() and plot_dendogram() functions, respectively
                                                                          # wt_cluster is the integer number referring to the Wild type cluster
    """Implementation of Hierarchical HotNet algorithm to find significantly altered subnetworks in each cluster"""

    # Setting environment
    start_time = time.time()
    if wt_cluster not in set(scoring_df['cluster']): raise ValueError('wt_cluster must be an integer number referring to an existing cluster')
    path = './results/hotnet_output'
    if not os.path.exists(f'{path}/data'): os.makedirs(f'{path}/data')
    if not os.path.exists(f'{path}/intermediate/network'): os.makedirs(f'{path}/intermediate/network')
    if not os.path.exists(f'{path}/results'): os.makedirs(f'{path}/results')
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
    lfc_df = lfc_df.transpose().reset_index(drop=False, names='display name')
    for cluster in cluster_dict.keys():
        for variant in cluster_dict[cluster]:
            if not os.path.exists(f'{path}/data/scores_{variant}.tsv'): pd.concat([lfc_df['display name'], lfc_df[variant].abs()], axis=1).to_csv(f'{path}/data/scores_{variant}.tsv', sep="\t", index=False, header=False)
            if not os.path.exists(f'{path}/intermediate/network_score_{variant}'): os.makedirs(f'{path}/intermediate/network_score_{variant}')

    # Constructing similarity matrices
    if os.path.exists(f'{path}/intermediate/network/similarity_matrix.h5') and os.path.exists(f'{path}/intermediate/network/beta.txt'): print('Similarity matrices already constructed')
    else:
        print('Constructing similarity matrices...')
        os.system(f'python ./hotnet/construct_similarity_matrix.py -i {path}/data/network_edges.tsv -o {path}/intermediate/network/similarity_matrix.h5 -bof {path}/intermediate/network/beta.txt')
    
    # Permuting networks
    cond = True
    for i in range(1, 9):
        if not os.path.exists(f'{path}/intermediate/network/edge_list_{i}.tsv'): cond = False
    if cond: print('Networks already permuted')
    else: 
        print('Permuting networks...')
        os.system(f'cp {path}/data/network_nodes.tsv {path}/intermediate/network/index_gene_0.tsv')
        os.system(f'cp {path}/data/network_edges.tsv {path}/intermediate/network/edge_list_0.tsv')
        for i in range(1, 5):  # Preserve connectivity of the observed graph
            os.system(f'python ./hotnet/permute_network.py -i {path}/intermediate/network/edge_list_0.tsv -s {i} -c -o {path}/intermediate/network/edge_list_{i}.tsv')
        for i in range(5, 9):  # Do not preserve connectivity of the observed graph
            os.system(f'python ./hotnet/permute_network.py -i {path}/intermediate/network/edge_list_0.tsv -s {i} -o {path}/intermediate/network/edge_list_{i}.tsv')

    # Permuting scores
    cond = True
    for cluster in cluster_dict.keys():
        for variant in cluster_dict[cluster]:
            if not os.path.exists(f'{path}/intermediate/network_score_{variant}/scores_0.tsv') or not os.path.exists(f'{path}/intermediate/network_score_{variant}/score_bins.tsv'): cond = False
            for i in range(1, num_perm + 1):
                if not os.path.exists(f'{path}/intermediate/network_score_{variant}/scores_{i}.tsv'): cond = False
    if cond: print('Scores already permuted')
    else: 
        print('Permuting scores...')
        for cluster in cluster_dict.keys():
            for variant in cluster_dict[cluster]:
                os.system(f'cp {path}/data/scores_{variant}.tsv {path}/intermediate/network_score_{variant}/scores_0.tsv')
                os.system(f'python ./hotnet/find_permutation_bins.py -gsf {path}/intermediate/network_score_{variant}/scores_0.tsv -igf {path}/data/network_nodes.tsv \
                -elf {path}/data/network_edges.tsv -ms  1000 -o {path}/intermediate/network_score_{variant}/score_bins.tsv')
                for i in range(1, num_perm + 1):
                    os.system(f'python ./hotnet/permute_scores.py -i  {path}/intermediate/network_score_{variant}/scores_0.tsv \
                    -bf {path}/intermediate/network_score_{variant}/score_bins.tsv -s  {i} -o {path}/intermediate/network_score_{variant}/scores_{i}.tsv')

    # Constructing hierarchies...
    cond = True
    for cluster in cluster_dict.keys():
        for variant in cluster_dict[cluster]:
            for i in range(num_perm + 1):
                if not os.path.exists(f'{path}/intermediate/network_score_{variant}/hierarchy_edge_list_{i}.tsv') or not os.path.exists(f'{path}/intermediate/network_score_{variant}/hierarchy_index_gene_{i}.tsv'): cond = False
    if cond: print('Hierarchies already constructed')
    else: 
        print('Constructing hierarchies...')
        for cluster in cluster_dict.keys():
            for variant in cluster_dict[cluster]:
                for i in range(num_perm + 1):
                    if not os.path.exists(f'{path}/intermediate/network_score_{variant}/hierarchy_edge_list_{i}.tsv') or not os.path.exists(f'{path}/intermediate/network_score_{variant}/hierarchy_index_gene_{i}.tsv'):
                        os.system(f'python ./hotnet/construct_hierarchy.py -smf  {path}/intermediate/network/similarity_matrix.h5 -igf  {path}/data/network_nodes.tsv \
                    -gsf  {path}/intermediate/network_score_{variant}/scores_{i}.tsv -helf {path}/intermediate/network_score_{variant}/hierarchy_edge_list_{i}.tsv -higf {path}/intermediate/network_score_{variant}/hierarchy_index_gene_{i}.tsv')

    # Processing hierarchies...
    cond = True
    for cluster in cluster_dict.keys():
        for variant in cluster_dict[cluster]:
            if not os.path.exists(f'{path}/results/clusters_network_scores_{variant}.tsv') or not os.path.exists(f'{path}/results/sizes_network_scores_{variant}.pdf'): cond = False
    if cond: print('Hierarchies already processed')
    else: 
        print('Processing hierarchies...')
        for cluster in cluster_dict.keys():
            for variant in cluster_dict[cluster]:
                os.system(f'python ./hotnet/process_hierarchies.py -oelf {path}/intermediate/network_score_{variant}/hierarchy_edge_list_0.tsv -oigf {path}/intermediate/network_score_{variant}/hierarchy_index_gene_0.tsv \
                -pelf $(for i in `seq {num_perm}`; do echo " {path}/intermediate/network_score_{variant}/hierarchy_edge_list_"$i".tsv "; done) -pigf $(for i in `seq {num_perm}`; do echo " {path}/intermediate/network_score_{variant}/hierarchy_index_gene_"$i".tsv "; done) \
                -lsb  1 -cf {path}/results/clusters_network_scores_{variant}.tsv -pl network scores_{variant} -pf {path}/results/sizes_network_scores_{variant}.pdf')

    # Performing consensus in each cluster...
    cond = True
    for cluster in cluster_dict.keys():
        if not os.path.exists(f'{path}/results/consensus_nodes_cluster_{cluster}.tsv') or not os.path.exists(f'{path}/results/consensus_edges_cluster_{cluster}.tsv'): cond = False
    if cond: print('Consensus already performed')
    else:
        print('Performing consensus...')
        for cluster in cluster_dict.keys():
            cluster_list = cluster_dict[cluster]
            cf = ' '.join([f'{path}/results/clusters_network_scores_{variant}.tsv' for variant in cluster_list])
            igf = ' '.join([f'{path}/data/network_nodes.tsv' for i in range(len(cluster_list))])
            elf = ' '.join([f'{path}/data/network_edges.tsv' for i in range(len(cluster_list))])
            n = ' '.join(['network' for i in range(len(cluster_list))])
            s = ' '.join([f'score_{variant}' for variant in cluster_list])
            os.system(f'python ./hotnet/perform_consensus.py -cf {cf} -igf {igf} -elf {elf} -n {n} -s {s} -t 2 -cnf {path}/results/consensus_nodes_cluster_{cluster}.tsv -cef {path}/results/consensus_edges_cluster_{cluster}.tsv')
            consensus_list = []
            consensus_df = pd.read_csv(f'{path}/results/consensus_nodes_cluster_{cluster}.tsv', sep='\t', header=None)
            for row in consensus_df.index:
                consensus_list += list(consensus_df.loc[row].dropna())
            pd.DataFrame({'consensus_nodes': consensus_list}).to_csv(f'{path}/results/consensus_nodes_cluster_{cluster}.txt', index=None)

    print(f'Consensus nodes and edges altered in each cluster saved in {path}/results/')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")