# Parameters

## scRNA-Seq Data Parameters
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| single cell data file* | Single-cell gene expression data file. Takes a [.h5ad file](https://github.com/scverse/anndata/issues/180#issuecomment-510451045) input. | No default value |
| cell type metadata field* | Field name in the sc/snRNA-Seq data containing cell type groupings. | No default value |

## Spatial Data Parameters
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| spatial data file* | Spatial gene expression data file. Takes a [.h5ad file](https://github.com/scverse/anndata/issues/180#issuecomment-510451045) input. | No default value |
| cluster metadata field | Field name in the spatial data containing cluster groupings. | No default value |

## "Genes to Plot" Parameter
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| genes to plot data file* | Input data file for genes to compare measured vs predicted data of. Reads the initial line of a [.gmt file](https://www.genepattern.org/file-formats-guide#GMT). | No default value |

## Deconvolution Data Parameters (OPTIONAL)
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| histological image file | A .zip file of a [zarr store](https://anndata.readthedocs.io/en/latest/generated/anndata.read_zarr.html) containing a histology image. <b><u>Important Note:</b></u> Only provide an input here to run the additional deconvolution step. | No default value |
| img metadata field | Name of the image layer in the <i>histological image file</i> object to process during deconvolution. | No default value |


## Output Filename Parameter
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| output file tag* | Output filename tag to be used as a suffix for this run's output. It is encouraged to change this parameter for every run you intend to do as to differentiate between outputs. The default uses the basenames of both the sc/snRNA-seq and spatial input datasets. | <single.cell.data.file_basename>_<spatial.data.file_basename> |

## Model Training Parameters
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| classification mode* | Mode of classification to carry out Tangram as either one of <b><u>Test-Train</b></u> or <b><u>Leave-One-Out Cross-Validation</b></u>. | Test-Train |
| training mode* | Mode of training for the Tangram algorithm as one of the following:  <u><b>Top N Genes</b></u> uses the top N differentially expressed genes stratified across cell types (you can set N in the next parameter, <i>num top DE genes</i>), <u><b>GMT File Input</u></b> uses an input [.gmt file](https://www.genepattern.org/file-formats-guide#GMT) where the first line's genes are used to train the model (you can input the file in <i>training genes input</i>), <u><b>All Genes</u></b> uses all genes shared between the two datasets. | Top N Genes |
| num top DE genes* | Integer value of top differentially expressed genes, stratified across cell types and shared between the two datasets, to use for training the prediction model. <u><b>Note:</b></u> This is only considered if <i>training mode</i> is set to <b>Top N Genes</b>. | 100 |
| training genes input | File input for the genes to be used in model training. Reads a single-line [.gmt file](https://www.genepattern.org/file-formats-guide#GMT). <b><u>Note:</b></u> This is only used if <i>training mode</i> is set to <b>GMT File Input</b>. | No default value |

## Algorithmic Parameters
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| alignment mode* | Mode of alignment or mapping. Use <u><b>cells</b></u> for mapping at a single-cell resolution, <b><u>clusters</b></u> for averaging single cells of the same cluster via the <i>cell type metadata field</i>, and <b><u>constrained</b></u> for deconvolution which adds a filter term to the loss function. | cells |
| alignment cell density* | Alignment cell density mode. Use <u><b>RNA count-based</b></u> if cell density is proportional to number of RNA molecules and <u><b>Uniform</b></u> if spatial voxels are at single cell resolution. | RNA count-based |
| num epochs* | Integer value for number of epochs, or iterations through the entire training dataset for the mapping process. | 1000 |

## Visualization Parameters
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| spot size* | Float value for diameter of gene plot spots. <u><b>Note:</b></u> This is only considered if there is no 'spatial' field in the spatial data input object (<i>spatial data file</i>). | 50 |
| scale factor* | Float value for gene plot scaling factor to map dots into space. <u><b>Note:</b></u> This is only considered if there is no 'spatial' field in the spatial data input object (<i>spatial data file</i>).| 0.1 |
| spatial opacity* | Float value for opacity in single-cell scanpy spatial plot of original data. | 0.7 |
| training opacity* | Float value of opacity for training score plots. | 0.5 |
| testing opacity* | Float value of opacity for testing score plots. | 0.5 |
| deconv opacity* | Float value for opacity in single-cell scanpy spatial plot of post-deconvolution cell plot. | 0.2 |
| deconv point size* | Float value for point size of single-cell scanpy spatial plot in post-deconvolution cell plot. | 0.4 |
| training bin num* | Integer value of bins in the training plot histogram. | 20 |
| testing bin num* | Integer value of bins in the testing plot histogram. | 20 |
| color mapping perc* | Percentile float (0.0 to 1.0, inclusive of both) value for setting color map range (a smaller value corresponds to the color map having a greater coverage of values). | 0.02 |

## Developer Parameters (Ignore)
<b>\* required field</b>

| Name | Description | Default Value |
---------|--------------|----------------|
| debug* | Developer debugging int (0 = no verbosity, > 0 for verbose). | 0 |
| verbose* | Developer verbosity int (0 = no verbosity, > 0 for verbose). | 0 |

