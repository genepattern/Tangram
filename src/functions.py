import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import tangram as tg
import anndata as ad

def h5ad_input(input):
    """
    Read in user input, a .h5ad file of either sc/snrna or spatial data.
    
    Args:
        input (str): filename of a .h5ad containing single cell or spatial data input

    Returns:
        AnnData object
    """

    # Removed .gct input code (available in previous commits)

    adata = sc.read_h5ad(input)
    return adata

def plot_initial(adata_sc, adata_sp, sc_celltype_field, sp_cluster_field, sc_alpha, umap_dotsize, sc_filename, umap_filename):
    """
    Outputs initial plots of sc/snrna and spatial AnnData objects to files.
    
    Args:
        adata_sc (AnnData): object of single-cell data
        adata_sp (AnnData): object of spatial data
        sc_celltype_field (str): name of .obs field in sc/snrna data with cell-type groupings
        sp_cluster_field (str): name of .obs field in spatial data with cluster groupings 
        sc_alpha (float): opacity value for single-cell data plot (lower is more transparent)
        umap_dotsize (int): umap dot size for spatial data plot
        sc_filename (str): output filename for single-cell data plot
        umap_filename (str): output filename for umap plot
    
    Returns:
        N/A
    """

    fig, ax = plt.subplots(figsize=(10, 5))
    sc.pl.spatial(adata_sp, color=sp_cluster_field, alpha=sc_alpha, frameon=False, show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(sc_filename)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    sc.pl.umap(adata_sc, color=sc_celltype_field, size=umap_dotsize, frameon=False, show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(umap_filename)
    plt.close(fig)

def gmt_input(filename):
    """
    Read in user input .gmt and convert to a list of genes as strings.
        
    Args:
        input (filename): filename of .gmt file with a single line (one set of genes)

    Returns:
        Genes as a list of strings
    """

    with open(filename, 'r') as file:
        line = file.readline()
        cells = line.strip().split('\t')
        genes = cells[2:]

    return genes

# TODO: Clarify in module and documentation the usage for cv and tt modes
def read_markers(adata_sc, mode, sc_celltype_field, n, gmt_file):
    """
    Only utilized in BOTH test-train and cv prediction modes.
    Evaluates training mode (use top n genes across cell types, .gmt file input, or all genes shared between the datasets) and returns training genes (all genes to be used in case of cv).
        
    Args:
        adata_sc (AnnData): object of single-cell data
        mode (str): test-train mode ("Top N Genes," "GMT File Input," or "All Genes")
        sc_celltype_field (str): name of .obs field in sc/snrna data with cell-type groupings
        n (int): for "Top N Genes," the value of N across cell types
        gmt_file (str): for "GMT File Input," the filename of the one-line .gmt

    Returns:
        List of training genes or None (in the case of "All Genes")
    """

    if (mode == "Top N Genes"):
        sc.tl.rank_genes_groups(adata_sc, groupby=sc_celltype_field, use_raw=False)
        markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:n, :]
        markers = list(np.unique(markers_df.melt().value.values))
    elif (mode == "GMT File Input"):
        markers = gmt_input(gmt_file)
    else:
       markers = None

    return markers

def preprocessing(adata_sc, adata_sp, markers):
    """
    Use Tangram's pp_adatas for preprocessing (only keeps shared genes, removes genes with 0-count in a dataset, etc).
    
    Args:
        adata_sc (AnnData): object of single cell data
        adata_sp (AnnData): object of spatial data
        markers (list of strs): training genes list as strings

    Returns:
        N/A (pp_adatas processes input objects)       
    """

    tg.pp_adatas(adata_sc, adata_sp, genes=markers)

def cv(adata_sc, adata_sp, mode, cv_mode, sc_celltype_field, device, num_epochs, verbose):
    """
    Performs cross-validation on the data.
    
    Args:
        adata_sc (AnnData): object of single cell data
        adata_sp (AnnData): object of spatial data
        mode (str): choice of either one of:
                "cells" (default, GPU recommended), 
                "clusters" (averages single cells belonging to same cluster),
                "constrained" (for deconvolution)
        cv_mode (str): either one of "loo" (leave-one-out) or "10fold"
        sc_celltype_field (str): for "constrained" mode, name of .obs field in sc/snrna data with cell-type groupings
        device (str): "cuda:0" (GPU) or "cpu"
        num_epochs (int): # of iterations for function mapping

    Returns:
        ad_ge_cv (AnnData): object of new spatial mapping
        df_test_genes (DataFrame): matrix of test genes' scores and sparsities 
    """

    cv_dict, ad_ge_cv, df_test_genes = tg.cross_val(adata_sc, adata_sp, 
                                     device=device, 
                                     mode=mode,
                                     cv_mode=cv_mode,
                                     num_epochs=num_epochs,
                                     cluster_label=sc_celltype_field,
                                     return_gene_pred=True,
                                     verbose=verbose)
    
    print("Average cross-validation scores:", cv_dict)

    return ad_ge_cv, df_test_genes

def find_alignment(adata_sc, adata_sp, mode, sc_celltype_field, density_prior, num_epochs, device):
    """
    Use Tangram's map_cells_to_space to create alignment map.

    Args:
        adata_sc (AnnData): object of single cell data
        adata_sp (AnnData): object of spatial data
        mode (str): choice of either one of:
                        "cells" (default, GPU recommended), 
                        "clusters" (averages single cells belonging to same cluster),
                        "constrained" (for deconvolution)
        sc_celltype_field (str): for "constrained" mode, name of .obs field in sc/snrna data with cell-type groupings
        density_prior (str): choice of either one of:
                                "uniform" (if spatial voxels at single-cell resolution),
                                "rna_count_based" (assumes cell density is proportional to # of RNA molecules),
                                for "constrained" mode, is calculated to be the fraction of cells per voxel        
        num_epochs (int): # of iterations for function mapping
        device (str): "cuda:0" (GPU) or "cpu"

    Returns:
        AnnData object containing mapping
        - cell-by-spot matrix X - probability of cell i to be in spot j
        - obs (dataframe) - metadata of single cells
        - var (dataframe) - metadata of spatial data
        - uns (dictionary) - dataframe with various info about training genes (saved as train_genes_df).
    """

    # Accounting for deconvolution
    if (mode == "constrained"):
        target_count = adata_sp.obs.cell_count.sum()
        density_prior = np.array(adata_sp.obs.cell_count) / target_count
    else:
        target_count = None

    # Accounting for cluster mode:
    if (mode != "clusters"):
        sc_celltype_field = None

    adata_map = tg.map_cells_to_space(adata_sc, adata_sp,
        mode=mode,
        target_count=target_count,
        density_prior=density_prior,
        num_epochs=num_epochs,
        cluster_label=sc_celltype_field,
        device=device
    )

    return adata_map

def map_scrna_to_space(adata_map, adata_sc, adata_sp, sc_celltype_field, celltype_map_filename, perc):
    """
    Use Tangram's project_cell_annotations to transfer annotation from mapping to space, then plot.
    
    Args:
        adata_map (AnnData): object of alignment map
        adata_sc (AnnData): object of single-cell data
        adata_sp (AnnData): object of spatial data
        sc_celltype_field (str): name of .obs field in sc/snrna data with cell-type groupings
        perc (float): colormap range
        celltype_map_filename (str): Filename to output the plot
        
    Returns:
        N/A
    """

    tg.project_cell_annotations(adata_map, adata_sp, annotation=sc_celltype_field)
    annotation_list = list(pd.unique(adata_sc.obs[sc_celltype_field]))
    tg.plot_cell_annotation_sc(adata_sp, annotation_list, perc)
    
    plt.savefig(celltype_map_filename)
    plt.close()

# Replacement function of tg.plot_training_scores() for saving output plot
def plot_training_replacement(adata_map, training_alpha, num_bins, filename):
    """
    Plots the 4-panel training diagnosis plot.

    Args:
        adata_map (AnnData): object containing alignment map
        num_bins (int): number of histogram bins
        alpha (float): plot opacity value (lower is more transparent)
        filename (str): name of file to output plots to

    Returns:
        N/A
    """

    fig, axs = plt.subplots(1, 4, figsize=(12, 3), sharey=True)
    df = adata_map.uns["train_genes_df"]
    axs_f = axs.flatten()

    axs_f[0].set_ylim([0.0, 1.0])
    for i in range(1, len(axs_f)):
        axs_f[i].set_xlim([0.0, 1.0])
        axs_f[i].set_ylim([0.0, 1.0])

    axs_f[0].set_title('Training scores for single genes')
    sns.histplot(data=df, y="train_score", bins=num_bins, ax=axs_f[0], color="coral")

    axs_f[1].set_title("score vs sparsity (single cells)")
    sns.scatterplot(
        data=df,
        y="train_score",
        x="sparsity_sc",
        ax=axs_f[1],
        alpha=training_alpha,
        color="coral",
    )

    axs_f[2].set_title("score vs sparsity (spatial)")
    sns.scatterplot(
        data=df,
        y="train_score",
        x="sparsity_sp",
        ax=axs_f[2],
        alpha=training_alpha,
        color="coral",
    )

    axs_f[3].set_title("score vs sparsity (sp - sc)")
    sns.scatterplot(
        data=df,
        y="train_score",
        x="sparsity_diff",
        ax=axs_f[3],
        alpha=training_alpha,
        color="coral",
    )

    plt.tight_layout()
    plt.savefig(filename)

# Replacement function of tg.plot_test_scores() for saving output plot
def plot_test_replacement(df_gene_score, bins, alpha, filename):
    """
    Plots gene level test scores with each gene's sparsity for mapping result.
    
    Args:
        df_gene_score (Pandas dataframe): returned by compare_spatial_geneexp(adata_ge, adata_sp, adata_sc); 
                       with "gene names" as the index and "score", "sparsity_sc", "sparsity_sp", "sparsity_diff" as the columns
        bins (int): number of histogram bins
        alpha (float): plot opacity value (lower is more transparent)
        filename (str): name of file to output plots to

    Returns:
        None
    """

    if not set(["score", "sparsity_sc", "sparsity_sp", "sparsity_diff"]).issubset(
        set(df_gene_score.columns)
    ):
        raise ValueError(
            "There are missing columns in df_gene_score. Run `compare_spatial_geneexp` with `adata_ge`, `adata_sp`, and `adata_sc` to produce complete dataframe input."
        )

    if "is_training" in df_gene_score.keys():
        df = df_gene_score[df_gene_score["is_training"] == False]
    else:
        df = df_gene_score

    df.rename({"score": "test_score"}, axis="columns", inplace=True)

    fig, axs = plt.subplots(1, 4, figsize=(12, 3), sharey=True)
    axs_f = axs.flatten()

    axs_f[0].set_ylim([0.0, 1.0])
    for i in range(1, len(axs_f)):
        axs_f[i].set_xlim([0.0, 1.0])
        axs_f[i].set_ylim([0.0, 1.0])

    sns.histplot(data=df, y="test_score", bins=bins, ax=axs_f[0])

    axs_f[1].set_title("score vs sparsity (single cells)")
    sns.scatterplot(data=df, y="test_score", x="sparsity_sc", ax=axs_f[1], alpha=alpha)

    axs_f[2].set_title("score vs sparsity (spatial)")
    sns.scatterplot(data=df, y="test_score", x="sparsity_sp", ax=axs_f[2], alpha=alpha)

    axs_f[3].set_title("score vs sparsity (sp - sc)")
    sns.scatterplot(
        data=df, y="test_score", x="sparsity_diff", ax=axs_f[3], alpha=alpha
    )
    plt.tight_layout()
    plt.savefig(filename)

def project(adata_map, adata_sc, predictions_filename):
    """
    Generates "new spatial data" using the now-mapped sc/snrna data.
        
    Args:
        adata_map (AnnData): object of spatial mapping (from find_alignment()) 
        adata_sc (AnnData): object of single-cell data
        predictions_filename (str): filename of new spatial data object to output

    Returns:
        AnnData object of new spatial mapping
    """

    ad_ge = tg.project_genes(adata_map, adata_sc)
    ad_ge.write(predictions_filename)

    return ad_ge

def compare_spatial(ad_ge, adata_sp, adata_sc, data_filename):
    """
    Compute similarity scores of all genes between measured and predicted values.
        
    Args:
        adata_sp (AnnData): object of original spatial data
        adata_sc (AnnData): object of single-cell data
        adata_ge (AnnData): object of new spatial data (from project())
        data_filename (str): filename of new spatial data object to output

    Returns:
        DataFrame of all genes' (train and test) similarity scores and sparsity statistics
    """
    df_all_genes = tg.compare_spatial_geneexp(ad_ge, adata_sp, adata_sc)
    df_all_genes.to_csv(data_filename, sep='\t')

    return df_all_genes

# Replacement function of tg.plot_auc() due to seaborn version issues
def plot_auc_replacement(df_all_genes, plot_filename, test_genes=None):
    """
      Plots auc curve which is used to evaluate model performance.
    
    Args:
        df_all_genes (DataFrame): matrix of all genes, train and test, and their sparsity statistics (from compare_spatial())
        plot_filename (str): filename to output auc plot to
        test_genes (list of strs): test genes, and if not given, test_genes will be set to genes where 'is_training' field is False

    Returns:
        N/A
    """
    
    metric_dict, ((pol_xs, pol_ys), (xs, ys)) = tg.eval_metric(df_all_genes, test_genes)
    
    fig = plt.figure()
    plt.figure(figsize=(6, 5))

    plt.plot(pol_xs, pol_ys, c='r')
    sns.scatterplot(data=pd.DataFrame({'x': xs, 'y': ys}), x="x", y="y", alpha=0.5)
        
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.gca().set_aspect(.5)
    plt.xlabel('score')
    plt.ylabel('spatial sparsity')
    plt.tick_params(axis='both', labelsize=8)
    plt.title('Prediction on test transcriptome')
    
    textstr = 'auc_score={}'.format(np.round(metric_dict['auc_score'], 3))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.3)

    plt.text(0.03, 0.1, textstr, fontsize=11, verticalalignment='top', bbox=props)
    plt.savefig(plot_filename)

def plot_genes(genes_filepath, adata_sp, ad_ge, plot_genes_outpath, perc):
    """
    Plots measured vs predicted transcriptomic gene profiles in space 
        
    Args:
        genes_filepath (str): filename of one-line GMT file to take in genes to plot
        adata_sp (AnnData): object of original spatial data
        adata_ge (AnnData): object of new spatial data (from project())
        plot_genes_outpath (str): filename to output gene predicted vs measured plots to
        perc (float): colormap range

    Returns:
        N/A
    """

    genes = gmt_input(genes_filepath)
    fig = tg.plot_genes_sc(genes, adata_measured=adata_sp, adata_predicted=ad_ge, perc=perc, return_figure=True)

    fig.savefig(plot_genes_outpath)

# TODO: Make parameter assignments more verbose
def execute_tangram_workflow(args):
    """
    Function to execute the Tangram workflow based on parsed arguments.

    Args:
        args (argparse.Namespace): parsed arguments from command line

    Returns:
        None
    """

    adata_sp = h5ad_input(args.sp)
    adata_sc = h5ad_input(args.sc)
    celltype_label = args.sc_celltype_field

    plot_initial(adata_sc, adata_sp, celltype_label, args.sp_cluster_field, args.spatial_alpha, args.umap_point_size, args.spatial_plot_filename, args.umap_plot_filename)

    markers = read_markers(adata_sc, args.training_mode, celltype_label, args.number_training_genes, args.marker_genes_input)

    preprocessing(adata_sc, adata_sp, markers)

    if (args.classification_mode == "Cross-Validation"):
        ad_ge, df_test_genes = cv(adata_sc, adata_sp, args.alignment_mode, args.cross_val_mode, args.device, args.num_epochs, args.verbose)
        plot_test_replacement(df_test_genes, args.test_bin_num, args.testing_alpha, args.test_plot_filename)
    else:
        adata_map = find_alignment(adata_sc, adata_sp, args.alignment_mode, celltype_label, args.density_prior, args.num_epochs, args.device)
        map_scrna_to_space(adata_map, adata_sc, adata_sp, celltype_label, args.celltype_plot_filename, args.perc)
        plot_training_replacement(adata_map, args.training_alpha, args.train_bin_num, args.training_plot_filename)
        ad_ge = project(adata_map, adata_sc, args.predictions_filename)

        df_all_genes = compare_spatial(ad_ge, adata_sp, adata_sc, args.similarity_scores_filename)
        plot_auc_replacement(df_all_genes, args.auc_plot_filename)

    plot_genes(args.genes_to_plot, adata_sp, ad_ge, args.genes_plot_filename, args.perc)
    
    return None
