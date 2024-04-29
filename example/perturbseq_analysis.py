import os 
import os 
import time
import numpy as np
import pandas as pd
import scanpy as sc
import pickle as pkl
import matplotlib.pyplot as plt
from gprofiler import GProfiler
from scipy.stats import ttest_ind
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from matplotlib.backends.backend_pdf import PdfPages
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from pydeseq2.default_inference import DefaultInference
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def filter_data(matrix, var_names, obs_names, mt_names):       # matrix is a dataframe whose head and index are genes and cells, respectively,
    """Filtering data; check comments below"""                 # var_names is a list containing the genes (HGNC nomenclature) assigned to each column,
    start_time = time.time()                                   # obs_names is a list containing the variants assigned to each cell (row),
    if not os.path.exists('./results'): os.mkdir('./results')  # and mt_names is a list containing the mitochondrial genes names

    # Checking parameters                                 
    if matrix.shape[0] != len(obs_names): raise ValueError('Number of rows in matrix does not match the length of variants list')
    if matrix.shape[1] != len(var_names): raise ValueError('Number of columns in matrix does not match the length of genes list')

    # Filtering matrix
    total_counts = list(matrix.sum(1))
    matrix = matrix.iloc[[index for index, count in enumerate(total_counts) if count > 7000], :]      # Filtering cells with <7000 counts
    obs_names = [variant for index, variant in enumerate(obs_names) if index in matrix.index]         # Updating variants list
    matrix = matrix.reset_index(drop=True)
    matrix.columns = var_names
    mt_genes_found = [mt_name for mt_name in mt_names if mt_name in list(var_names)]                  # Finding mitochondrial genes in matrix
    mt_counts = list(matrix[mt_genes_found].sum(1))
    pct_mt_counts = [(i / j) * 100 for i, j in zip(mt_counts, total_counts)]
    matrix = matrix.iloc[[index for index, count in enumerate(pct_mt_counts) if count < 20], :]       # Filtering cells with >20% mitochondrial counts
    total_counts = list(matrix.sum(1))                                                                # Updating total counts, variants list, and % mito counts
    obs_names, pct_mt_counts = zip(*[(tuple[0], tuple[1]) for index, tuple in enumerate(zip(obs_names, pct_mt_counts)) if index in matrix.index])
    matrix = matrix.reset_index(drop=True)

    # Creating AnnData
    adata = sc.AnnData(matrix)
    adata.var['var_names'] = list(var_names)
    adata.obs['obs_names'] = list(obs_names)
    adata.obs['total_counts'] = list(total_counts)
    adata.obs['pct_mt_counts'] = list(pct_mt_counts)
    adata.obs['n_genes_counts'] = list(matrix.astype(bool).sum(axis=1))

    # Filtering AnnData
    sc.pp.filter_cells(adata, min_genes=200)                # Filtering cells with <200 genes
    sc.pp.filter_genes(adata, min_cells=3)                  # Filtering genes present in <3 cells
    sc.pl.highest_expr_genes(adata, n_top=20, save='.png')  # Plotting 20 most highly variable genes
    sc.pp.downsample_counts(adata, counts_per_cell=50000)   # Downsampling cells with >50,000 counts to <=50,000 counts

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
    sc.pp.neighbors(adata, n_neighbors=5)
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
    sc.pp.neighbors(adata, n_neighbors=5)
    sc.tl.louvain(adata, key_added='louvain')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['louvain'], s=5, save = '.png')
    os.system('mv ./figures/umap.png ./results/louvain_genes.png')
    os.system('rm -r ./figures')

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
    matrix = pd.DataFrame(adata.X)
    matrix.columns = adata.var['var_names']
    # adata.obs.index = range(len(adata.obs.index))
    matrix.index = adata.obs.index
    if not os.path.exists('./results/diff_analysis'): os.mkdir('./results/diff_analysis')
    
    for variant in set(adata.obs['obs_names']):
        if variant != 'WT' and variant != 'Null':
            # Selecting variants
            matrix_temp = matrix[(adata.obs['obs_names'] == 'WT') | (adata.obs['obs_names'] == variant)].reset_index(drop=True)
            obs_temp = adata.obs[(adata.obs['obs_names'] == 'WT') | (adata.obs['obs_names'] == variant)].reset_index(drop=True)
        
            # Initializing read counts modeling
            inference = DefaultInference(n_cpus=8)
            dds = DeseqDataSet(counts=matrix_temp, metadata=obs_temp, design_factors='obs_names', refit_cooks=True, inference=inference)
            dds.deseq2()                                     # Fitting dispersions and log fold changes
        
            # Statistical analysis
            stat_res = DeseqStats(dds, inference=inference)  # Compute p-values using Wald tests and adjusted p-values for differential expresion
            stat_res.summary()                               # Running the whole statistical analysis, cooks filtering, and multiple testing adjustement
            result_dict[variant] = stat_res
    
    with open('./results/diff_analysis/diff_analysis.pkl', 'wb') as file:
        pkl.dump(result_dict, file)
    print('Display results of each variant vs WT doing returned_dict[variant].results_df')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")
    return result_dict  # result_dict is a dict: keys are variants and values are DEA results with regards to WT:

def plot_dea(adata_dict):  # result_dict is the dict obtained from the diff_analysis function
    """Plotting the DEA"""
    start_time = time.time()
    for variant in adata_dict.keys():
        stat_res = adata_dict[variant]
        stat_res.plot_MA(s=5)
       #plt.gcf().set_size_inches(8,6)
        plt.gcf().set_dpi(100)
        plt.title(variant + '_WT')
        plt.tight_layout()
        plt.savefig('./results/diff_analysis/' + variant + '_WT.png', dpi=150)
    print('For a higher resolution of a particular variant, type returned_dict[variant].plot_MA(s=5)')
    print("Execution time:", round(time.time() - start_time, 3), "seconds")