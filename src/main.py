## python script to run main function 
## Add any argument parsing, function calling, here

import argparse
from functions import * 

### main function to parse arguments
## calls from functions.py to run functions
def main():

  # File input arguments

  # Spatial data file input (.gct):
  parser.add_argument("--sp", "--spatial", help="spatial data file input"
                      + " Valid file format(s): .gct")
  # Single-cell data file input (.gct):
  parser.add_argument("--sc", "--scrna", help="single-cell data file input"
                      + " Valid file format(s): .gct")


  # Tangram parameter arguments
  # TODO: remove nargs and const, keeping them now for likely smoother debugging - feel free to change  

  # Classification mode ("test-train" or "loocv"):
  parser.add_argument("-m", "--classification_mode", help="tangram classification mode (test-train or loocv)",
                      nargs="?", const=1, default="test-train", choices=["test-train","loocv"])
  
  # UMAP point size argument, used for one plot (initial umap); integer value:
  parser.add_argument("--umap_point_size", help="umap point size in single-cell initial scanpy plot",
                    nargs="?", const=1, type=int, default=10)
  
  # Alpha value for scanpy spatial plot opacity; float value ranging from 0.0 to 1.0, inclusive:
  parser.add_argument("--spatial_alpha", help="alpha value in single-cell scanpy spatial plots",
                    nargs="?", const=1, type=float, default=0.7)
  
  # Alpha value for training score plot opacity; float value ranging from 0.0 to 1.0, inclusive:
  parser.add_argument("--training_alpha", help="alpha value in training score plots",
                    nargs="?", const=1, type=float, default=0.5)

  # Boolean value for whether to use top n differentially expressed genes as training genes:
  parser.add_argument("--use_top_n", help="boolean for whether to use top n differentially expressed" 
                    + "genes for training", nargs="?", const=1, type=bool, default=True)
  
  # Number of top differentially expressed genes to use as training genes, only considered if use_top_n is true:
  parser.add_argument("-n", "--number_training_genes", help="number of top differentially expressed genes to"
                    + "use for training", nargs="?", const=1, type=int, default=100)
  
  # File input for manual training marker gene selection for preprocessing step, only considered if use_top_n is false:
  # TODO: decide on necessity of parameter and file format (simple tsv of names?)
  parser.add_argument("--marker_genes_input", help="file input in case of manual training marker gene selection",
                    nargs="?", const=1)

  # Alignment mode argument, accepted values are either "cluster", "cell", or "constrained"; gpu pref if "cell":
  parser.add_argument("--alignment_mode", help="tangram alignment mode (cluster, cell, or constrained)",
                      nargs="?", const=1, default="cluster", choices=["cluster","cell","constrained"])
 
  # Alignment cell density argument, accepted values are either "rna_count_based" (cell density prop. to number of RNA molecules) 
  # or "uniform" (if spatial voxels at single cell resolution):
  parser.add_argument("--density_prior", help="tangram alignment cell density within each spatial voxel (uniform or rna_count_based)",
                      nargs="?", const=1, default="rna_count_based", choices=["rna_count_based","uniform"])
  
  # Color-mapping percent for plots argument, accepted values are floats between 0 and 1:
  # TODO: look into how to potentially obtain a somewhat optimal value of this for users, might be useful to investigate
  # some of the commented-out "robust" argument code in tg; could potentially calculate this value for users based on outliers
  parser.add_argument("--perc", help="tangram range of color-mapping for plots (adjust to remove outliers)",
                      nargs="?", const=1, default=0.02, type=float)

  # TODO: Make educated decision about default # of epochs, current default value set same as tutorial
  parser.add_argument("--num_epochs", help="Number of iterations for function mapping", type=int, default=500)
  
  # TODO: Potentially adjust this so that choices are "cpu" and "cuda", for user friendliness
  parser.add_argument("--device", help="Device to use (cpu or cuda:0 for GPU)", default="cpu", choices=["cpu", "cuda:0"])
  
  # File input of genes to plot measured vs predicted for:
  # TODO: decide on file format (simple tsv of names?)
  parser.add_argument("--genes_to_plot", help="file input of genes to plot measured vs predicted for",
                    nargs="?", const=1)
  # Annotation used when projecting and plotting cell annotations in deconvolution step
  parser.add_argument("--annotation", help="Annotation parameter for spatial mapping", default="cell_subclass", type = "str")


  # Developer & verbosity arguments
  # Module verbosity argument, either True (1) or False (0), False by default
  parser.add_argument("-v", "--verbose", help="module verbosity flag",
                      nargs="?", const=1, default=0, type=int)
  # Module debug argument, either True (1) or False (0), False by default
  parser.add_argument("-d", "--debug", help="module debug flag",
                      nargs="?", const=1, default=0, type=int)
  
  return None




## Parsing from command line, and running the script. 
if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  args = parser.parse_args()

  if (args.debug):
    print("Debugging on.\n")

