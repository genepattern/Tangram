## python script to include required functions for running
## add any required functions here
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import tangram as tg
import anndata as ad

def genes_list_from_file(genes_input):
    """
    Read a tab-separated file of gene names and convert it into a list of strings.

    Args:
        file_path (str): Path to the tab-separated file containing gene names.

    Returns:
        list: A list of gene names as strings.
    """
    # Read gene names from the tab-separated file
    with open(genes_input, 'r') as file:
        genes_list = [line.strip() for line in file]
    return genes_list

  
def plot_auc_replacement(df_all_genes, test_genes=None):
    """
      A replacement function of the tangram-sc library's plot_auc.
      Plots auc curve which is used to evaluate model performance.
    
    Args:
        df_all_genes (Pandas dataframe): returned by compare_spatial_geneexp(adata_ge, adata_sp); 
        test_genes (list): list of test genes, if not given, test_genes will be set to genes where 'is_training' field is False

    Returns:
        None
    """
    metric_dict, ((pol_xs, pol_ys), (xs, ys)) = tg.eval_metric(df_all_genes, test_genes)
    
    fig = plt.figure()
    plt.figure(figsize=(6, 5))

    plt.plot(pol_xs, pol_ys, c='r')
    sns.scatterplot(data=pd.DataFrame({'x': xs, 'y': ys}), alpha=0.5)
        
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
    plt.text(0.03, 0.1, textstr, fontsize=11, verticalalignment='top', bbox=props);


def read_data(sc_input, spatial_input):
  """
    Read in user inputs, scrnaseq data and spatial data.
    
    Args:
        sc_input (filename): GCT file containing single cell data input
        spatial_data (filename): GCT file containing spatial data input

    Returns:
        Dataframes for respective input files (sc_data, spatial_data)
  """
  sc_data = pd.read_csv(sc_input, sep='\t', skiprows=2)
  spatial_data = pd.read_csv(spatial_input, sep='\t', skiprows=2)
  return sc_data, spatial_data


def convert_gct_to_adata(gct_df):
    """
    Convert GCT (Gene Cluster Text) input dataframe to AnnData
    
    Args:
        gct_file

    Returns:
        AnnData object
    """
    gene_names = gct_df.columns[2:]
    sample_names = gct_df.iloc[:, 0]
    expression_data = gct_df.iloc[:, 2:].T # remove gene and sample names
  
    adata = ad.AnnData(X=expression_data.values, obs=sample_names, var=gene_names)
    return adata

# TODO: adjust to alternatively use use_top_n and number_training_genes, prob create helper function for this
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
  # return None



def find_alignments(adata_sc, adata_sp, mode, density, num_epochs, device, cluster_label=None):
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
  adata_map = tg.map_cells_to_space(adata_sc, adata_sp, mode, cluster_label, 
                                    density, num_epochs, device)

  return adata_map


def map_scrna_to_space(adata_map, adata_sc, adata_sp, annotation, perc=0.02):
  """
    Use Tangram's project_cell_annotations to transfer annotation from mapping to space, then plot
    
    Args:
        adata_map (AnnData): alignment map
        adata_sc (AnnData): single-cell data
        adata_sp (AnnData): spatial data
        annotation (str)
        perc(double): colormap range
        

    Returns:
        Plots with spatial mapping
        
  """
  
  tg.project_cell_annotations(adata_map, adata_sp, annotation)
  # Should plotting be in this function or separate instead?
  annotation_list = list(pd.unique(adata_sc.obs[annotation])) # based on tutorial code
  tg.plot_cell_annotation_sc(adata_sp, annotation_list,perc)
  # return None


def plot_genes_sc(genes, adata_sp, adata_ge, perc=0.02):
  """
    Use Tangram's plot_genes_sc to create gene expression maps from the two spatial AnnData on the genes
    
    Args:
        genes (list of strs): genes list
        adata_sp (AnnData): measured spatial data
        adata_ge (AnnData): predicted spatial data
        perc(double): colormap range
        

    Returns:
        Plots with spatial mapping
        
  """
  tg.plot_genes_sc(genes, perc, adata_measured=adata_sp, adata_predicted=adata_ge)
  

def execute_tangram_workflow(args):
    """
    Function to execute the Tangram workflow based on parsed arguments.

    Args:
        args (argparse.Namespace): Parsed arguments from command line

    Returns:
        None
    """
    # Read data
    sc_data, spatial_data = read_data(args.scrna, args.spatial)

    # Convert both single-cell and spatial data to AnnData objects
    adata_sc = convert_gct_to_adata(sc_data)
    adata_sp = convert_gct_to_adata(spatial_data)

    # Preprocessing
    markers = genes_list_from_file(args.marker_genes_input)
    preprocessing(adata_sc, adata_sp, markers)

    # Find alignments
    adata_map = find_alignments(adata_sc, adata_sp, args.alignment_mode, args.density_prior,
                                args.num_epochs, args.device)

    # Note: currently the same color-mapping percent is assumed between the cell type maps and following actual/predicted cell gene pattern plot 
    # Map single-cell data to spatial space
    map_scrna_to_space(adata_map, adata_sc, adata_sp, args.annotation, args.perc)

    # Plot genes
    if args.genes_to_plot:
        genes_to_plot = genes_list_from_file(args.genes_to_plot)
        plot_genes_sc(genes_to_plot, adata_sp, adata_map, args.perc)


    return None
