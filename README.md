# Tangram
#### Omar Halawa (ohalawa@ucsd.edu) & Julia Kononova (jkononova@ucsd.edu) of the GenePattern Team @ Mesirov Lab - UCSD

## Summary
The following repository is a GPU-enabled GenePattern module of [Tangram](https://github.com/broadinstitute/Tangram) (see [publication](https://www.nature.com/articles/s41592-021-01264-7)), a deep learning method for mapping sc/snRNA-seq data to various forms of spatial data collected from the same anatomical region or tissue type.   

   ![Tangram](other/tangram_description.png)

It can be used for various purposes, including but not limited to: 
- extending gene throughput ![Cell Types](other/cell_types.png)
- obtaining a **spatial localization of cell types**
- correcting low-quality data
- performing single-cell **deconvolution**  <img src="/other/deconv.png" alt="deconv" width="1000"/>


## Implementation   

It was written in Python 3 and uses a Singularity image built from its [Docker counterpart](https://hub.docker.com/layers/omarhalawa/tangram/v1.3/images/sha256-77a432f4d5ebe023e0ed073cb945600dd2f316e14b20aa27cc1f4077cf88f840?context=repo). 

Documentation on all module parameters can be found [here](/docs/parameters.md).
All source files, including input and output datasets for replicating runs from the original tutorial notebooks ([tutorial_tangram_with_squidpy.ipynb](https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_with_squidpy.ipynb) & [tutorial_tangram_without_squidpy.ipynb](https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_without_squidpy.ipynb)), can be found [here](/data/).   

## References   

Biancalani, Tommaso, et al. **“Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram.”** _Nature Methods_, vol. 18, no. 11, 28 Oct. 2021, pp. 1352–1362.</p>
