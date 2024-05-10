import argparse
from functions import * 

def main():

  # File input arguments
  # Spatial data file input (.gct/.h5ad):
  parser.add_argument("--sp", "--spatial", help="spatial data file input"
                      + " Valid file format(s): .h5ad")
  # Single-cell data file input (.gct/.h5ad):
  parser.add_argument("--sc", "--scrna", help="single-cell data file input"
                      + " Valid file format(s): .h5ad")
  # GMT file input of genes to plot measured vs predicted for (.gmt):
  parser.add_argument("--genes_to_plot", help="file input of genes to plot measured vs predicted for"
                      + " Valid file format(s): .gmt")

  # Output filename arguments
  # Filename for UMAP plot of initial single-cell data
  parser.add_argument("--umap_plot_filename", help="filename for initial umap plot output")
  # Filename for plot of initial spatial data 
  parser.add_argument("--spatial_plot_filename", help="filename for initial spatial plot output")
  # Filename for training plot
  parser.add_argument("--train_plot_filename", help="filename for training plot")
  # Filename for testing plot
  parser.add_argument("--test_plot_filename", help="filename for testing plot")
  # Filename for plotted genes
  parser.add_argument("--genes_plot_filename", help="filename for plotted genes")
  # Filename for auc plot
  parser.add_argument("--auc_plot_filename", help="filename for auc plot")
  # Filename for celltype plot
  parser.add_argument("--celltype_plot_filename", help="filename for celltype plot")
  # Filename for anndata object of predictions
  parser.add_argument("--predictions_filename", help="filename for anndata predictions output")
  # Filename for dataframe object of similarity scores of genes
  parser.add_argument("--similarity_scores_filename", help="filename for output dataframe of similarity scores")

  # Tangram parameter arguments
  # Classification mode ("Test-Train" or "Cross-Validation"):
  parser.add_argument("--classification_mode", help="tangram classification mode (test-train or xval)",
                      nargs="?", const=1, default="Test-Train", choices=["Test-Train","Cross-Validation"])
  # Cross-Validation mode:
  parser.add_argument("--cross_val_mode", help="CV mode, either 'loo' for loocv or '10fold' as of the tg 1.0.4",
                      nargs="?", const=1, default="Test-Train", choices=["loo","10fold"])
  
  # Single-cell cell-type field name:
  parser.add_argument("--sc_celltype_field", help="name of .obs field in sc/snrna data for the cell-type groupings",
                    type=str, default="cell_subclass")
  # Spatial cluster field name:
  parser.add_argument("--sp_cluster_field", help="name of .obs field in spatial data for the cluster groupings",
                    type=str, default="clusters")
  
  # UMAP point size argument, used for one plot (initial umap); integer value:
  parser.add_argument("--umap_point_size", help="umap point size in single-cell initial scanpy plot",
                    type=int, default=10)
  
  # Alpha value for scanpy spatial plot opacity; float value ranging from 0.0 to 1.0, inclusive:
  parser.add_argument("--spatial_alpha", help="alpha value in single-cell scanpy spatial plots",
                    type=float, default=0.7)
  
  # Boolean value for whether to use top n differentially expressed genes as training genes:
  parser.add_argument("--training_mode", help="choice for whether to use top n differentially expressed" 
                    + "genes shared for training, gmt gene set, or all shared genes", nargs="?", const=1, default="Top N Genes", 
                    choices=["Top N Genes", "GMT File Input", "All Genes"])

  # Alpha value for training score plot opacity; float value ranging from 0.0 to 1.0, inclusive:
  parser.add_argument("--training_alpha", help="alpha value in training score plots",
                    type=float, default=0.5)
  # Alpha value for testing score plot opacity; float value ranging from 0.0 to 1.0, inclusive:
  parser.add_argument("--testing_alpha", help="alpha value in testing score plots",
                    type=float, default=0.5)

  # Number of top differentially expressed genes to use as training genes, only considered if use_top_n is true:
  parser.add_argument("-n", "--number_training_genes", help="number of top differentially expressed genes to"
                    + "use for training", type=int, default=100)
  
  # File input for manual training marker gene selection, only considered if use_top_n is false:
  parser.add_argument("--marker_genes_input", help="file input in case of manual training marker gene selection")

  parser.add_argument("--train_bin_num", help="number of bins in the training histogram plot",
                    type=int, default=20)
  parser.add_argument("--test_bin_num", help="number of bins in the testing histogram plot",
                    type=int, default=20)

  # Alignment mode argument, accepted values are: "clusters", "cells", or "constrained"; gpu pref if "cell":
  parser.add_argument("--alignment_mode", help="tangram alignment mode (cluster, cell, or constrained)",
                      default="cluster", choices=["clusters","cells","constrained"])
  
  # TODO: explain that this parameter is only utilized in the cases of no deconvolution (not "constrained" mode)
  # Alignment cell density argument, accepted values are either "rna_count_based" (cell density prop. to number of RNA molecules) 
  # or "uniform" (if spatial voxels at single cell resolution):
  parser.add_argument("--density_prior", help="tangram alignment cell density within each spatial voxel (uniform or rna_count_based)",
                      default="rna_count_based", choices=["rna_count_based","uniform"])
    
  # Color-mapping percent for plots argument, accepted values are floats between 0 and 1:
  # TODO: look into how to potentially obtain a somewhat optimal value of this for users, might be useful to investigate
  # some of the commented-out "robust" argument code in tg; could potentially calculate this value for users based on outliers
  parser.add_argument("--perc", help="tangram range of color-mapping for plots (adjust to remove outliers)",
                      default=0.02, type=float)

  # TODO: Make educated decision about default # of epochs, current default value set same as tutorial
  parser.add_argument("--num_epochs", help="Number of iterations for function mapping", type=int, default=500)

  # TODO: Potentially adjust this so that choices are "cpu" and "cuda", for user friendliness
  parser.add_argument("--device", help="Device to use (cpu or cuda:0 for GPU)", default="cpu", choices=["cpu", "cuda:0"])

  # Developer & verbosity arguments
  # Module verbosity argument, either True (1) or False (0), False by default
  parser.add_argument("-v", "--verbose", help="module verbosity flag", default=0, type=int)
  # Module debug argument, either True (1) or False (0), False by default
  parser.add_argument("-d", "--debug", help="module debug flag", default=0, type=int)
  
  return None

# TODO: Handle errors more gracefully (filepath not found, etc..)
# Parsing from command line and running the script
if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  main()
  args = parser.parse_args()

  if (args.debug):
    print("Debugging on.\n")
  
  execute_tangram_workflow(args)
