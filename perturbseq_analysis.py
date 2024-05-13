import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import pickle as pkl
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
    sc.pl.highest_expr_genes(adata, n_top=20, save='.png')
    sc.pp.highly_variable_genes(adata, n_top_genes=int(len(adata.var.index) * 0.2), flavor='seurat_v3', inplace=True)
    adata = adata[:,adata.var['highly_variable']]
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
        
            # Initializing read counts modeling
            dds = DeseqDataSet(counts=pd.DataFrame(adata_temp.X, columns = adata.var['var_names']), metadata=adata_temp.obs, design_factors='obs_names', refit_cooks=True, inference=inference)
            dds.deseq2()                                     # Fitting dispersions and log fold changes
        
            # Statistical analysis
            stat_res = DeseqStats(dds, inference=inference)  # Compute p-values using Wald tests and adjusted p-values for differential expresion
            stat_res.summary()                               # Running the whole statistical analysis, cooks filtering, and multiple testing adjustement
            result_dict[variant] = stat_res
    
    with open('./results/diff_analysis/diff_analysis.pkl', 'wb') as file:
        pkl.dump(result_dict, file)
    print('Display results of each variant vs WT doing returned_dict[variant].results_df')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_dict  # result_dict is a dict: keys are variants and values are DEA results with regards to WT

def plot_dea(adata_dict):  # result_dict is the dict obtained from the diff_analysis function
    """Plotting the DEA, and saving info only about LFC and p values"""
    start_time = time.time()
    result_df = pd.DataFrame()
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
        result_df = pd.concat([result_df, stat_res.results_df[['log2FoldChange', 'padj']].rename(columns={'log2FoldChange': 'LFC_' + variant, 'padj': 'padj_' + variant})], axis=1)

    result_df.fillna(1, inplace=True)  # Set NaN values in p value columns to 1
    result_df.to_csv('./results/diff_analysis/diff_analysis.csv')
    print('For a higher resolution plot of a particular variant, type returned_dict[variant].plot_MA(s=5)')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return filtered_dict  # filtered_dict is a dict: keys are variants and values are differential expressed genes with their LFC values