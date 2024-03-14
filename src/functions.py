## python script to include required functions for running
## add any required functions here
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import tangram as tg

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


def read_data():
  ## add any parameters as necessary
  ## reads user inputs, scrnaseq data and spatial data.
  return None


def convert_gct_to_adata():
  ## add any parameters as necessary
  ## converts gct file to adata format
  return None


def preprocessing():
  ## add any parameters as necessary
  return None



def find_alignments():
  ## add any parameters as necessary
  return None


def map_scrna_to_space():
  ## add any parameters as necessary
  return None


def plot_genes_sc(genes):
  ## add any parameters as necessary
  ## takes a list of genes and calls tangram plot_genes_sc to generate
  ## spatial 
  return None

