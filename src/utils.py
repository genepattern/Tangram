import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import skimage
import scanpy as sc
import squidpy as sq
import tangram as tg
import anndata as ad
import zipfile

def h5ad_input(filepath):
    """
    Read in user input, a .h5ad file of either sc/snrna or spatial data.
    
    Args:
        filepath (str): filename of a .h5ad containing single cell or spatial data input

    Returns:
        AnnData object
    """

    # Removed .gct input code (available in previous commits)

    adata = sc.read_h5ad(filepath)
    return adata

def img_input(filepath):
    """
    Read in user input, a zipped file of a zarr store of a squidpy ImageContainer object.
    
    Args:
        filepath (str): filename of a zarr store zip

    Returns:
        Squidpy ImageContainer object (or None if no filepath is provided)
    """

    if (filepath == None or filepath == ""):
        return None

    with zipfile.ZipFile(filepath, 'r') as zip_ref:
        zip_ref.extractall("/".join(filepath.split("/")[:-1]))
    
    # This always works as we guarantee a .zip input in the module
    img = sq.im.ImageContainer.load(filepath[:-4])

    return img

def gmt_input(filepath):
    """
    Read in user input .gmt and convert to a list of genes as strings.
        
    Args:
        filepath (str): filename of .gmt file with a single line (one set of genes)

    Returns:
        Genes as a list of strings
    """

    with open(filepath, 'r') as file:
        line = file.readline()
        cells = line.strip().split('\t')
        genes = cells[2:]

    return genes

# Removed due to bugs caused by inconsistency in data inputs (without_squidpy.ipynb cannot run)
# def plot_initial_umap(adata_sc, sc_celltype_field,  umap_dotsize, umap_filename):
#     """
#     Outputs initial umap plot of sc/snRNA AnnData object.
    
#     Args:
#         adata_sc (AnnData): object of single-cell data
#         sc_celltype_field (str): name of .obs field in sc/snrna data with cell type groupings
#         umap_dotsize (int): umap dot size for spatial data plot
#         umap_filename (str): output filename for umap plot
    
#     Returns:
#         N/A
#     """

#     fig, ax = plt.subplots(figsize=(10, 5))
#     sc.pl.umap(adata_sc, color=sc_celltype_field, size=umap_dotsize, frameon=False, show=False, ax=ax)
#     plt.tight_layout()
#     plt.savefig(umap_filename)
#     plt.close(fig)

def plot_initial_spatial(adata_sp, sp_cluster_field, sc_alpha, sc_filename):
    """
    Outputs initial plot of spatial AnnData object.
    
    Args:
        adata_sp (AnnData): object of spatial data
        sp_cluster_field (str): name of .obs field in spatial data with cluster groupings 
        sc_alpha (float): opacity value for single-cell data plot (lower is more transparent)
        sc_filename (str): output filename for single-cell data plot
    
    Returns:
        N/A
    """

    fig, ax = plt.subplots(figsize=(10, 5))
    sc.pl.spatial(adata_sp, color=sp_cluster_field, alpha=sc_alpha, frameon=False, show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(sc_filename)
    plt.close(fig)

def read_markers(adata_sc, mode, sc_celltype_field, n, gmt_file):
    """
    Utilized in BOTH test-train and cv prediction modes.
    Evaluates training mode (use top n genes across cell types, .gmt file input, or all genes shared between the datasets) and returns training genes (all genes to be used in case of cv).
        
    Args:
        adata_sc (AnnData): object of single-cell data
        mode (str): test-train mode ("top_n," "gmt_input," or "all_genes")
        sc_celltype_field (str): name of .obs field in sc/snrna data with cell type groupings
        n (int): for "top_n," the value of N across cell types
        gmt_file (str): for "gmt_input," the filename of the one-line .gmt

    Returns:
        List of training genes or None (in the case of "All Genes")
    """

    if (mode == "top_n"):
        sc.tl.rank_genes_groups(adata_sc, groupby=sc_celltype_field, use_raw=False)
        markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:n, :]
        markers = list(np.unique(markers_df.melt().value.values))
    elif (mode == "gmt_input"):
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

def cv(adata_sc, adata_sp, mode, cv_mode, sc_celltype_field, device, num_epochs, verbose, data_filename, matrix_filename):
    """
    Performs cross-validation on the data.

    Args:
        adata_sc (AnnData): Single cell data object.
        adata_sp (AnnData): Spatial data object.
        mode (str): Mode of operation - "cells", "clusters", or "constrained".
        cv_mode (str): Cross-validation mode - "loo" (leave-one-out) or "10fold".
        sc_celltype_field (str): Field in .obs for cell type groupings (for "constrained" mode).
        device (str): "cuda:0" (GPU) or "cpu".
        num_epochs (int): Number of iterations for function mapping.
        verbose (bool): Verbosity of the output.
        data_filename (str): Filename for saving new spatial mapping (.h5ad).
        matrix_filename (str): Filename for saving test score matrix (.csv).

    Returns:
        ad_ge_cv (AnnData): New spatial mapping (only if "loo" mode).
        df_test_genes (DataFrame): Matrix of test gene scores and sparsities.
    """

    common_params = {
        'adata_sc': adata_sc,
        'adata_sp': adata_sp,
        'device': device,
        'mode': mode,
        'cv_mode': cv_mode,
        'num_epochs': num_epochs,
        'cluster_label': sc_celltype_field,
        'return_gene_pred': True,
        'verbose': verbose
    }

    if cv_mode == "10fold":
        cv_dict, df_test_genes = tg.cross_val(**common_params)
        print("Average cross-validation scores:", cv_dict)
        return df_test_genes
    
    if cv_mode == "loo":
        cv_dict, ad_ge_cv, df_test_genes = tg.cross_val(**common_params)
        print("Average cross-validation scores:", cv_dict)
        
        ad_ge_cv.write(data_filename)
        ad_ge_cv.var.sort_values(by='test_score', ascending=False).to_csv(matrix_filename)
        
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
        sc_celltype_field (str): for "constrained" mode, name of .obs field in sc/snrna data with cell type groupings
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

def map_scrna_to_space(adata_map, adata_sc, adata_sp, sc_celltype_field, celltype_map_filename, perc, spot_size, scale_factor):
    """
    Use Tangram's project_cell_annotations to transfer annotation from mapping to space, then plot.
    
    Args:
        adata_map (AnnData): object of alignment map
        adata_sc (AnnData): object of single-cell data
        adata_sp (AnnData): object of spatial data
        sc_celltype_field (str): name of .obs field in sc/snrna data with cell type groupings
        perc (float): colormap range
        celltype_map_filename (str): Filename to output the plot
        spot_size (int): diameter of plot spots; only used if no "spatial" field in adata_sp
        scale_factor (float): scaling factor used to map from coordinate space to pixel space; only used if no "spatial" field in adata_sp

    Returns:
        N/A
    """

    if 'spatial' in adata_sp.uns.keys():
        spot_size=None
        scale_factor=None

    tg.project_cell_annotations(adata_map, adata_sp, annotation=sc_celltype_field)
    annotation_list = list(pd.unique(adata_sc.obs[sc_celltype_field]))
    tg.plot_cell_annotation_sc(adata_sp, annotation_list, perc=perc, spot_size=spot_size, scale_factor=scale_factor)
    
    plt.tight_layout()
    plt.savefig(celltype_map_filename)
    plt.close()

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
    df_all_genes.to_csv(data_filename)

    return df_all_genes

def segment_cells(img, layer):
    """
    Use squidpy to segment cells on the corresponding histology image. This is done to know how many cells are present in each voxel.
    
    Args:
        img (ImageContainer): squidpy ImageContainer object of histology data
        layer (str): name of field containing the image layer to be processed

    Returns:
        N/A
    """
    sq.im.process(img=img, layer=layer, method="smooth")
    sq.im.segment(
        img=img,
        layer="image_smooth",
        method="watershed",
        channel=0,
    )

    # Commented out as unable to properly segment with the following variables without prior information of image from user
    # inset_y = 1500
    # inset_x = 1700
    # inset_sy = 400
    # inset_sx = 500

    # num_figs = 3 if sp_cluster_field != None else 2
    # fig, axs = plt.subplots(1, num_figs, figsize=(30, 10))

    # axs[0].imshow(
    #     img[layer][inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx, 0, 0]
    #     / 65536,
    #     interpolation="none",
    # )
    # axs[0].grid(False)
    # axs[0].set_xticks([])
    # axs[0].set_yticks([])
    # axs[0].set_title("DAPI", fontdict={"fontsize": 20})

    # crop = img["segmented_watershed"][
    #     inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx
    # ].values.squeeze(-1)
    # crop = skimage.segmentation.relabel_sequential(crop)[0]
    # cmap = plt.cm.plasma
    # cmap.set_under(color="black")
    # axs[1].imshow(crop, interpolation="none", cmap=cmap, vmin=0.001)
    # axs[1].grid(False)
    # axs[1].set_xticks([])
    # axs[1].set_yticks([])
    # axs[1].set_title("Nucleous segmentation", fontdict={"fontsize": 20})

    # if sp_cluster_field != None and sp_cluster_field != "":
    #     sc.pl.spatial(
    #         adata_sp, color=sp_cluster_field, alpha=0.7, frameon=False, show=False, ax=axs[0], title=""
    #     )
    #     axs[2].set_title("Clusters", fontdict={"fontsize": 20})
    #     axs[2].axes.xaxis.label.set_visible(False)
    #     axs[2].axes.yaxis.label.set_visible(False)

def extract_segmentation_features(adata_sp, img, layer, plot_filename):
    """
    Use Tangram's project_cell_annotations to transfer annotation from mapping to space, then plot.
    
    Args:
        adata_sp (AnnData): object of spatial data
        img (ImageContainer): squidpy ImageContainer object of histology data
        layer (str): name of field containing the image layer to be processed
        plot_filename (str): filepath to output cell density plot to
        
    Returns:
        N/A
    """

    # define image layer to use for segmentation
    features_kwargs = {
        "segmentation": {
            "label_layer": "segmented_watershed",
            "props": ["label", "centroid"],
            "channels": [1, 2],
        }
    }
    # calculate segmentation features
    sq.im.calculate_image_features(
        adata_sp,
        img,
        layer=layer,
        key_added="image_features",
        features_kwargs=features_kwargs,
        features="segmentation",
        mask_circle=True,
    )

    adata_sp.obs["cell_count"] = adata_sp.obsm["image_features"]["segmentation_label"]

    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.spatial(adata_sp, color="cell_count", frameon=False, show=False, ax=ax)
    plt.tight_layout()
    plt.savefig(plot_filename)
    plt.close(fig)

def deconvolve(adata_sc, adata_sp, adata_map, sc_celltype_field, deconv_alpha, deconv_spot_size,
               plot_filename, deconv_segmentation_data_filename, deconv_cell_type_counts_data_filename, 
               deconv_cell_type_annotations_data_filename):
    """
    Create cell segmentation of spatial data and output it, and then, perform deconvolution and plot cell types in space.
    
    Args:
        adata_sc (AnnData): object of single-cell data
        adata_sp (AnnData): object of spatial data
        adata_map (AnnData): object containing alignment map
        sc_celltype_field (str): name of .obs field in sc/snrna data with cell type groupings
        deconv_alpha (float): opacity value in deconvolution spatial plot
        deconv_spot_size (float):  spot size in deconvolution spatial plot
        plot_filename (str): filepath to output cell type plot following deconvolution
        deconv_segmentation_data_filename (str): csv filepath to output segmentation results
        deconv_cell_type_counts_data_filename (str): csv filepath to output cell type count results
        deconv_cell_annotations_data_filename (str): csv filepath to output cell annotation results
        
    Returns:
        N/A
    """
    
    tg.create_segment_cell_df(adata_sp)
    adata_sp.uns["tangram_cell_segmentation"].to_csv(deconv_segmentation_data_filename, index=False)

    tg.count_cell_annotations(adata_map, adata_sc, adata_sp, annotation=sc_celltype_field)
    adata_sp.obsm["tangram_ct_count"].to_csv(deconv_cell_type_counts_data_filename)

    adata_segment = tg.deconvolve_cell_annotations(adata_sp)
    adata_segment.obs.to_csv(deconv_cell_type_annotations_data_filename, index=False)

    fig, ax = plt.subplots(1, 1, figsize=(30, 30))
    sc.pl.spatial(
        adata_segment,
        color="cluster",
        size=deconv_spot_size,
        show=False,
        frameon=False,
        alpha_img=deconv_alpha,
        legend_fontsize=20,
        ax=ax,
    )

    plt.tight_layout()
    plt.savefig(plot_filename)

# Replacement function of tg.plot_training_scores() for saving output plot, also outputs csv of training data statistics
def plot_training_replacement(adata_map, training_alpha, num_bins, plot_filename, data_filename):
    """
    Plots the 4-panel training diagnosis plot.

    Args:
        adata_map (AnnData): object containing alignment map
        num_bins (int): number of histogram bins
        alpha (float): plot opacity value (lower is more transparent)
        plot_filename (str): name of file to output plots to
        data_filename (str): name of file to output training data matrix to

    Returns:
        N/A
    """

    df = adata_map.uns["train_genes_df"]
    df.to_csv(data_filename)

    fig, axs = plt.subplots(1, 4, figsize=(12, 3), sharey=True)
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
    plt.savefig(plot_filename)

# Replacement function of tg.plot_test_scores() for saving output plot
def plot_test_replacement(df_gene_score, bins, alpha, plot_filename):
    """
    Plots gene level test scores with each gene's sparsity for mapping result.
    
    Args:
        df_gene_score (Pandas dataframe): returned by compare_spatial_geneexp(adata_ge, adata_sp, adata_sc); 
                       with "gene names" as the index and "score", "sparsity_sc", "sparsity_sp", "sparsity_diff" as the columns
        bins (int): number of histogram bins
        alpha (float): plot opacity value (lower is more transparent)
        plot_filename (str): name of file to output plots to

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
    plt.savefig(plot_filename)

def project(adata_map, adata_sc, predictions_filename):
    """
    Generates "new spatial data" using the now-mapped sc/snrna data and outputs it as .h5ad.
        
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
    plt.tight_layout()
    plt.savefig(plot_filename)

def plot_genes(genes_filepath, adata_sp, ad_ge, plot_genes_outpath, perc, spot_size, scale_factor):
    """
    Plots measured vs predicted transcriptomic gene profiles in space 
        
    Args:
        genes_filepath (str): filename of one-line GMT file to take in genes to plot
        adata_sp (AnnData): object of original spatial data
        adata_ge (AnnData): object of new spatial data (from project())
        plot_genes_outpath (str): filename to output gene predicted vs measured plots to
        perc (float): colormap range
        spot_size (int): diameter of plot spots; only used if no "spatial" field in adata_sp
        scale_factor (float): scaling factor used to map from coordinate space to pixel space; only used if no "spatial" field in adata_sp

    Returns:
        N/A
    """

    if 'spatial' in adata_sp.uns.keys():
        spot_size=None
        scale_factor=None

    genes = gmt_input(genes_filepath)
    fig = tg.plot_genes_sc(genes, adata_measured=adata_sp, adata_predicted=ad_ge, perc=perc, return_figure=True, spot_size=spot_size, scale_factor=scale_factor)

    fig.savefig(plot_genes_outpath)

# TODO: Maybe make parameter assignments more verbose
def execute_tangram_workflow(args):
    """
    Function to execute the Tangram workflow based on parsed arguments.

    Args:
        args (argparse.Namespace): parsed arguments from command line

    Returns:
        None
    """

    adata_sc = h5ad_input(args.sc)
    adata_sp = h5ad_input(args.sp)

    if (args.sp_cluster_field != None and args.sp_cluster_field != ""):
        plot_initial_spatial(adata_sp, args.sp_cluster_field, args.spatial_alpha, f"initial_spatial_plot_{args.master_tag}.png")

    # Returns None if no image is provided
    img = img_input(args.img)
    if (args.verbose and img==None):
        print("No histological image provided.")

    markers = read_markers(adata_sc, args.training_mode, args.sc_celltype_field, args.num_training_genes, args.marker_genes_input)
    preprocessing(adata_sc, adata_sp, markers)

    if (args.classification_mode == "Cross-Validation"):
        cv_results = cv(adata_sc, adata_sp, args.alignment_mode, args.cross_val_mode, args.sc_celltype_field, args.device, args.num_epochs, 
                                  args.verbose, f"cv_predictions_{args.master_tag}.h5ad", f"testing_data_{args.master_tag}.csv")

        # Accounting for the different outputs with different cross-validation modes
        if args.cross_val_mode == "10fold":
            return
        else:
            plot_test_replacement(cv_results[-1], args.test_bin_num, args.testing_alpha, f"testing_plot_{args.master_tag}.png")
            ad_ge = cv_results[0]

    else:
        if (img):
            # Performing cell segmentation if an image is provided
            segment_cells(img, args.img_layer)
            extract_segmentation_features(adata_sp, img, args.img_layer, f"deconv_cell_density_plot_{args.master_tag}.png")

        adata_map = find_alignment(adata_sc, adata_sp, args.alignment_mode, args.sc_celltype_field, args.density_prior, args.num_epochs, args.device)
        map_scrna_to_space(adata_map, adata_sc, adata_sp, args.sc_celltype_field, f"cell_type_plot_{args.master_tag}.png", args.perc, args.spot_size, args.scale_factor)
        plot_training_replacement(adata_map, args.training_alpha, args.train_bin_num, f"training_plot_{args.master_tag}.png", f"training_data_{args.master_tag}.csv")
        ad_ge = project(adata_map, adata_sc, f"predictions_{args.master_tag}.h5ad")

        if (img):
            # Performing deconvolution if the input image is provided
            deconvolve(adata_sc, adata_sp, adata_map, args.sc_celltype_field, args.deconv_alpha, args.deconv_spot_size, 
                       f"deconv_cell_type_plot_{args.master_tag}.png", f"deconv_segmentation_{args.master_tag}.csv", 
                       f"deconv_cell_type_counts_{args.master_tag}.csv", f"deconv_cell_type_annotations_{args.master_tag}.csv")

        df_all_genes = compare_spatial(ad_ge, adata_sp, adata_sc, f"similarity_scores_{args.master_tag}.csv")
        plot_auc_replacement(df_all_genes, f"auc_plot_{args.master_tag}.png")

    plot_genes(args.genes_to_plot, adata_sp, ad_ge, f"gene_plots_{args.master_tag}.png", args.perc, args.spot_size, args.scale_factor)
    return None