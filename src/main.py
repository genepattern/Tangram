## python script to run main function 
## Add any argument parsing, function calling, here

import argparse
from functions import * 

### main function to parse arguments
## calls from functions.py to run functions
def main():

  # Adding arguments to script for Tangram file inputs, options, & debugging
  parser = ap.ArgumentParser(description='Tangram')

  # File input arguments
  # Spatial data file input (.gct):
  parser.add_argument("--sp", "--spatial", help="spatial data file input"
                      + " Valid file format(s): .gct")
  # Single-cell data file input (.gct):
  parser.add_argument("--sc", "--scrna", help="single-cell data file input"
                      + " Valid file format(s): .gct")

  # TODO: remove nargs and const, keeping them now for likely smoother debugging - feel free to change
  # Tangram parameter arguments
  # Alignment mode argument, accepted values are either "cluster", "cell", or "constrained"; gpu pref if "cell":
  # TODOC: see notebook for more details on these to have a concise parameter description in module
  parser.add_argument("-m", "--mode", help="tangram alignment mode (cluster, cell, or constrained)",
                      nargs="?", const=1, default="cluster", choices=["cluster","cell","constrained"])

  # TODO: add for training data what is likely three argparse arguments: one as a bool argument for whether or not to use
  # the top n differentially expressed genes, one as the n-value for top differentially expressed genes to use, and an argument
  # in case the user has a specific set of genes they want to use as training data and opted for F in the first arg (investigate data type for this)
 
  # Alignment cell density argument, accepted values are either "rna_count_based" (cell density prop. to number of RNA molecules) 
  # or "uniform" (if spatial voxels at single cell resolution):
  # TODO: investigate the np.array passed-in usage of this in the with-squidpy notebook (deconvolving part), likely not a MVP issue though
  parser.add_argument("--density_prior", help="tangram alignment cell density within each spatial voxel (uniform or rna_count_based)",
                      nargs="?", const=1, default="rna_count_based", choices=["rna_count_based","uniform"])
  
  # Color-mapping percent for plots argument, accepted values are floats between 0 and 1:
  # REF: change name of this param to something better? Essentially a substitute for the "perc" args set to 0.02 in the notebooks
  # TODO: look into how to potentially obtain a somewhat optimal value of this for users, might be useful to investigate
  # some of the commented-out "robust" argument code in tg; could potentially calculate this value for users based on outliers
  parser.add_argument("--plot_percent", help="tangram range of color-mapping for plots (adjust to remove outliers)",
                      nargs="?", const=1, default="cluster", type=float)

  # TODO: add an argument for a list of genes to plot individually using plot_genes_sc(), see cell 15-17 in the with-squidpy notebook
  # Unsure of whether to use string input param or file input - need to investigate

  
  # Developer & verbosity arguments
  # RFE: keep as ints (0/1 for F/T to be handled by module dropdown as "True"/"False") or other method? 
  # Module verbosity argument, either True or False, False by default
  parser.add_argument("-v", "--verbose", help="module verbosity flag",
                      nargs="?", const=1, default=False, type=int)
  # Module debug argument, either True or False, False by default
  parser.add_argument("-d", "--debug", help="module debug flag",
                      nargs="?", const=1, default=False, type=int)

  args = parser.parse_args()

  if (args.debug):
    print("Debugging on.\n")
  
  ## add any parameters as necessary
  return None




## Parsing from command line, and running the script. 
if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  args = parser.parse_args()
