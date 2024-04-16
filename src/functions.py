## python script to include required functions for running
## add any required functions here
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import tangram as tg
import anndata as ad

def read_data(input):
  """
    Read in user input, either a .gct or .h5ad file of scrna/spatial data.
    
    Args:
        input (filename): .gct/.h5ad file containing single cell or spatial data input

    Returns:
        AnnData object
  """
  
  if (input.endswith(".gct")):
    gct_df = pd.read_csv(input, sep='\t', skiprows=2)
    gene_names = gct_df.columns[2:]
    sample_names = gct_df.iloc[:, 0]
    expression_data = gct_df.iloc[:, 2:].T # remove gene and sample names
    adata = ad.AnnData(X=expression_data.values, obs=sample_names, var=gene_names)

  elif (input.endswith(".h5ad")):
    adata = sc.read_h5ad(input)

  return adata

def plot_initial(adata_sc, adata_sp, sc_alpha, umap_dotsize, sc_filename, umap_filename):
    """
    Outputs initial plots of sc and spatial anndata and saves them to files
    
    Args:
        adata_sc: anndata object of single-cell data
        adata_sp: anndata object of spatial data
        sc_alpha: alpha value for single-cell data plot 
        umap_dotsize: umap dot size for spatial data plot
        sc_filename: filename for single-cell data plot
        umap_filename: filename for umap plot
    
    Returns:
        N/A
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    sc.pl.spatial(adata_sp, color="cluster", alpha=sc_alpha, frameon=False, show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(sc_filename)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    sc.pl.umap(adata_sc, color="cell_subclass", size=umap_dotsize, frameon=False, show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(umap_filename)
    plt.close(fig)
    
def gmt_input(filename):
    """
        Read in user input .gmt and convert to a list of genes as strings
        
        Args:
            input (filename): .gmt file with a single line (one gene set)

        Returns:
            Genes as a list of strings
    """
    with open(filename, 'r') as file:
        line = file.readline()
        cells = line.strip().split('\t')
        genes = cells[2:]

    return genes

def read_markers(adata_sc, mode, n, gmt_file):
    """
        Evaluates mode (use top n genes across cell types or gmt file input) and returns training genes
        
        Args:
            adata_sc: single-cell anndata object
            mode: "Use Top N Genes" or "Use GMT File Input"
            n: For "Use Top N Genes", the value of N across cell types
            gmt_file: For "Use GMT File Input", the path to the one-line gmt

        Returns:
            List of training genes
    """

    if (mode == "Use Top N Genes"):
        sc.tl.rank_genes_groups(adata_sc, groupby="cell_subclass", use_raw=False)
        markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:n, :]
        markers = list(np.unique(markers_df.melt().value.values))
    else:
        markers = gmt_input(gmt_file)

    return markers

def preprocessing(adata_sc, adata_st, markers):
  """
    Use Tangram's pp_adatas for preprocessing
    
    Args:
        adata_sc (AnnData): single cell data
        adata_sp (AnnData): spatial data
        markers (list of strs): training genes list

    Returns:
        AnnData object containing mapping
        - cell-by-spot matrix X - probability of cell i to be in spot j
        - obs (dataframe) - metadata of single cells
        - var (dataframe) - metadata of spatial data
        - uns (dictionary) - dataframe with various info about training genes (saved as train_genes_df).
        
  """
  tg.pp_adatas(adata_sc, adata_st, genes=markers)

def find_alignment(adata_sc_in, adata_sp_in, mode_in, density_in, num_epochs_in, device_in, cluster_label_in=None):
  """
    Use Tangram's map_cells_to_space to create alignment map
    
    Args:
        adata_sc (AnnData): single cell data
        adata_sp (AnnData): spatial data
        mode (str): "cells" (default, GPU recommended) or "clusters" (averages single cells belonging to same cluster)
        cluster_label (str): annotations, for cluster-mode mapping only
        density_prior (str): "uniform" (if spatial voxels at single-cell resolution)
                            or "rna_count_based" (assumes cell density is proportional to # RNA molecules)          
        num_epochs (int): # of iterations for function mapping
        device (str): "cuda:0" (GPU) or "cpu"

    Returns:
        AnnData object containing mapping
        - cell-by-spot matrix X - probability of cell i to be in spot j
        - obs (dataframe) - metadata of single cells
        - var (dataframe) - metadata of spatial data
        - uns (dictionary) - dataframe with various info about training genes (saved as train_genes_df).
        
  """

  adata_map = tg.map_cells_to_space(adata_sc_in, adata_sp_in,
    mode=mode_in,
    density_prior=density_in,
    num_epochs=num_epochs_in,
    device=device_in,
    cluster_label = cluster_label_in
  )

  return adata_map

def map_scrna_to_space(adata_map, adata_sc, adata_sp, annotation_type, celltype_map_filename, perc):
    """
    Use Tangram's project_cell_annotations to transfer annotation from mapping to space, then plot
    
    Args:
        adata_map (AnnData): alignment map
        adata_sc (AnnData): single-cell data
        adata_sp (AnnData): spatial data
        annotation (str)
        perc (double): colormap range
        celltype_map_filename (str): Filename to save the plot (optional)
        
    Returns:
        Plots with spatial mapping
        
    """

    tg.project_cell_annotations(adata_map, adata_sp, annotation=annotation_type)
    annotation_list = list(pd.unique(adata_sc.obs[annotation_type]))
    tg.plot_cell_annotation_sc(adata_sp, annotation_list,perc)
    
    # Save the plot
    plt.savefig(celltype_map_filename)
    plt.close()

    return annotation_list

def plot_training_replacement(adata_map, training_alpha, num_bins, filename):
  """
  Plots the 4-panel training diagnosis plot

  Args:
      adata_map (AnnData):
      bins (int or string): Optional. Default is 10.
      alpha (float): Optional. Ranges from 0-1, and controls the opacity. Default is 0.7.

  Returns:
      None
  """
  fig, axs = plt.subplots(1, 4, figsize=(12, 3), sharey=True)
  df = adata_map.uns["train_genes_df"]
  axs_f = axs.flatten()

  # set limits for axis
  axs_f[0].set_ylim([0.0, 1.0])
  for i in range(1, len(axs_f)):
      axs_f[i].set_xlim([0.0, 1.0])
      axs_f[i].set_ylim([0.0, 1.0])

  #     axs_f[0].set_title('Training scores for single genes')
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

def project(adata_map, adata_sc, predictions_filename):
  ad_ge = tg.project_genes(adata_map, adata_sc)
  ad_ge.write(predictions_filename)
  return ad_ge

def plot_auc_replacement(ad_ge, adata_st, adata_sc, filename, test_genes=None):
    """
      A replacement function of the tangram-sc library's plot_auc.
      Plots auc curve which is used to evaluate model performance.
    
    Args:
        ad_ge (AnnData): anndata object of predicted spatial data
        adata_st: anndata object of input spatial data
        adata_sc: anndata object of input scrna data
        test_genes (list): list of test genes, if not given, test_genes will be set to genes where 'is_training' field is False

    Returns:
        None
    """

    df_all_genes = tg.compare_spatial_geneexp(ad_ge, adata_st, adata_sc)
    
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
    # place a text box in upper left in axes coords
    plt.text(0.03, 0.1, textstr, fontsize=11, verticalalignment='top', bbox=props)
    plt.savefig(filename)

def plot_genes(genes_filepath, adata_sp, ad_ge, plot_genes_outpath, perc=0.02):
  genes = gmt_input(genes_filepath)
  fig = tg.plot_genes_sc(genes, adata_measured=adata_sp, adata_predicted=ad_ge, perc=perc, return_figure=True)
  fig.savefig(plot_genes_outpath)

def execute_tangram_workflow(args):
    """
    Function to execute the Tangram workflow based on parsed arguments.

    Args:
        args (argparse.Namespace): Parsed arguments from command line

    Returns:
        None
    """
    ## Parsing from command line, and running the script. 

    adata_sp = read_data(args.sp)
    adata_sc = read_data(args.sc)

    plot_initial(adata_sc,adata_sp, args.spatial_alpha, args.umap_point_size, args.spatial_plot_filename, args.umap_plot_filename)

    markers = read_markers(adata_sc, args.training_mode, args.number_training_genes, args.marker_genes_input)

    preprocessing(adata_sc, adata_sp, markers)

    adata_map = find_alignment(adata_sc, adata_sp, mode_in=args.alignment_mode, density_in=args.density_prior, num_epochs_in=args.num_epochs, device_in=args.device, cluster_label_in=None)

    map_scrna_to_space(adata_map, adata_sc, adata_sp, args.annotation_type, args.celltype_plot_filename, args.perc)

    # TODO; handle the histogram bins input more smoothly
    plot_training_replacement(adata_map, args.training_alpha, 20, args.training_plot_filename)

    ad_ge = project(adata_map, adata_sc, args.predictions_filename)

    plot_genes(args.genes_to_plot, adata_sp, ad_ge, args.genes_plot_filename, args.perc)

    plot_auc_replacement(ad_ge, adata_sp, adata_sc, args.auc_plot_filename)
    
    return None
