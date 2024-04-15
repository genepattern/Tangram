import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
import tangram as tg

sc.logging.print_header()
print(f"squidpy=={sq.__version__}")


if __name__ == "__main__":
    adata_st = sq.datasets.visium_fluo_adata_crop()
    adata_st = adata_st[
        adata_st.obs.cluster.isin([f"Cortex_{i}" for i in np.arange(1, 5)])
    ].copy()
    img = sq.datasets.visium_fluo_image_crop()

    adata_sc = sq.datasets.sc_mouse_cortex()
    adata_st.obs

    fig, axs = plt.subplots(1, 2, figsize=(20, 5))
    sc.pl.spatial(
        adata_st, color="cluster", alpha=0.7, frameon=False, show=False, ax=axs[0]
    )
    sc.pl.umap(
        adata_sc, color="cell_subclass", size=10, frameon=False, show=False, ax=axs[1]
    )
    plt.tight_layout()

    sc.tl.rank_genes_groups(adata_sc, groupby="cell_subclass", use_raw=False)
    markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
    markers = list(np.unique(markers_df.melt().value.values))
    len(markers)
    
    tg.pp_adatas(adata_sc, adata_st, genes=markers)

    ad_map = tg.map_cells_to_space(adata_sc, adata_st,
    mode="cells",
#     mode="clusters",
#     cluster_label='cell_subclass',  # .obs field w cell types
    density_prior='rna_count_based',
    num_epochs=500,
    # device="cuda:0",
    device='cpu',
    )
    tg.project_cell_annotations(ad_map, adata_st, annotation="cell_subclass")
    annotation_list = list(pd.unique(adata_sc.obs['cell_subclass']))
    tg.plot_cell_annotation_sc(adata_st, annotation_list,perc=0.02)
    ad_map.uns['train_genes_df']
    ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
    genes = ['rragb', 'trim17', 'eno1b']
    ad_map.uns['train_genes_df'].loc[genes]
    tg.plot_genes_sc(genes, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)
    genes=['loc102633833', 'gm5700', 'gm8292']
    tg.plot_genes_sc(genes, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)
    (ad_ge.var.is_training == False).sum()
    df_all_genes = tg.compare_spatial_geneexp(ad_ge, adata_st, adata_sc)
    tg.plot_auc(df_all_genes)
    genes=['tfap2b', 'zic4']
    tg.plot_genes_sc(genes, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)
    genes = ['cd34', 'rasal1']
    tg.plot_genes_sc(genes, adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)
    genes = ['gm33027', 'gm5431']
    tg.plot_genes_sc(genes[:5], adata_measured=adata_st, adata_predicted=ad_ge, perc=0.02)
    sq.im.process(img=img, layer="image", method="smooth")
    sq.im.segment(
        img=img,
        layer="image_smooth",
        method="watershed",
        channel=0,
    )

    inset_y = 1500
    inset_x = 1700
    inset_sy = 400
    inset_sx = 500

    fig, axs = plt.subplots(1, 3, figsize=(30, 10))
    sc.pl.spatial(
        adata_st, color="cluster", alpha=0.7, frameon=False, show=False, ax=axs[0], title=""
    )
    axs[0].set_title("Clusters", fontdict={"fontsize": 20})
    sf = adata_st.uns["spatial"]["V1_Adult_Mouse_Brain_Coronal_Section_2"]["scalefactors"][
        "tissue_hires_scalef"
    ]
    rect = mpl.patches.Rectangle(
        (inset_y * sf, inset_x * sf),
        width=inset_sx * sf,
        height=inset_sy * sf,
        ec="yellow",
        lw=4,
        fill=False,
    )
    axs[0].add_patch(rect)

    axs[0].axes.xaxis.label.set_visible(False)
    axs[0].axes.yaxis.label.set_visible(False)

    axs[1].imshow(
        img["image"][inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx, 0, 0]
        / 65536,
        interpolation="none",
    )
    axs[1].grid(False)
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    axs[1].set_title("DAPI", fontdict={"fontsize": 20})

    crop = img["segmented_watershed"][
        inset_y : inset_y + inset_sy, inset_x : inset_x + inset_sx
    ].values.squeeze(-1)
    crop = skimage.segmentation.relabel_sequential(crop)[0]
    cmap = plt.cm.plasma
    cmap.set_under(color="black")
    axs[2].imshow(crop, interpolation="none", cmap=cmap, vmin=0.001)
    axs[2].grid(False)
    axs[2].set_xticks([])
    axs[2].set_yticks([])
    axs[2].set_title("Nucleous segmentation", fontdict={"fontsize": 20})
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
        adata_st,
        img,
        layer="image",
        key_added="image_features",
        features_kwargs=features_kwargs,
        features="segmentation",
        mask_circle=True,
    )
    adata_st.obs["cell_count"] = adata_st.obsm["image_features"]["segmentation_label"]
    sc.pl.spatial(adata_st, color=["cluster", "cell_count"], frameon=False)
    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_st,
        mode="constrained",
        target_count=adata_st.obs.cell_count.sum(),
        density_prior=np.array(adata_st.obs.cell_count) / adata_st.obs.cell_count.sum(),
        num_epochs=1000,
    #     device="cuda:0",
        device='cpu',
    )
    tg.project_cell_annotations(ad_map, adata_st, annotation="cell_subclass")
    annotation_list = list(pd.unique(adata_sc.obs['cell_subclass']))
    tg.plot_cell_annotation_sc(adata_st, annotation_list, perc=0.02)
    ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
    df_all_genes = tg.compare_spatial_geneexp(ad_ge, adata_st, adata_sc)
    tg.plot_auc(df_all_genes)
    tg.create_segment_cell_df(adata_st)
    adata_st.uns["tangram_cell_segmentation"].head()
    tg.count_cell_annotations(
    ad_map,
    adata_sc,
    adata_st,
    annotation="cell_subclass",
    )
    adata_st.obsm["tangram_ct_count"].head()
    adata_segment = tg.deconvolve_cell_annotations(adata_st)
    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    sc.pl.spatial(
        adata_segment,
        color="cluster",
        size=0.4,
        show=False,
        frameon=False,
        alpha_img=0.2,
        legend_fontsize=20,
        ax=ax,
    )
    
