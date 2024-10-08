# -*- coding: utf-8 -*-
"""Figure 4 Immune Organization/iPT/PT/iTAL/TAL Score HK2841.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1bN8FT-jZXxQE_dTgv_zudg46d2sMS2rt
"""

!pip install --quiet scanpy
!pip install --quiet leidenalg
!pip install --quiet squidpy

import csv
import anndata as ad
import gzip
import os
import scipy.io
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from google.colab import drive
import leidenalg as la
from pathlib import Path
import squidpy as sq
import scipy

from google.colab import drive
drive.mount('/content/drive')

adata3= sc.read_h5ad('/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/all_Cosmx_cleaned_neighbored.h5ad')

adata3.obs['iPT_subcluster_ME_20um'].value_counts()

sc.pl.umap(adata3, color= 'iPT_subcluster_ME')

adata3_subset=adata3[adata3.obs['iPT_subcluster_ME']!='outlier']
sc.pl.umap(adata3_subset, color = 'iPT_subcluster_ME')

adata3=adata3[adata3.obs["sample"]=="HK3035"]

adata3.obs['fov']=adata3.obs['fov'].astype(int)
#adata3.obs['fov']=adata3.obs['fov'].astype(str)

##this is for iPT_subcluster 4
fov=[74,75,76,77,78,79,80,81,
     82,83,84,85,86,87,88,89,
     90,91,92,93,94,95,96,97,
     98,99,100,101,102,103,104,105,
     106,107,108,109,110,111,112,113]
mask=adata3.obs['fov'].isin(fov)
adata3=adata3[mask]
adata3

##this is for iPT/iTAL intersection in the HK3035
fov=[6,13,20,28,36,44, 52, 60, 68,76,84]
mask=adata3.obs['fov'].isin(fov)
adata3=adata3[mask]
adata3

import numpy as np

# Create the new column 'for_plotting' based on the condition
adata3.obs['for_plotting'] = np.where(adata3.obs['iPT_subcluster_ME_20um'] == 'iPT ME 4',
                                      'iPT ME 4', 'outlier')
print(adata3.obs['for_plotting'].value_counts())

adata3.obs['for_plotting_iTAL'] = np.where(adata3.obs['iTAL_subcluster_ME_20um'] == 'iTAL ME 1',
                                      'iTAL ME 1', 'outlier')
print(adata3.obs['for_plotting_iTAL'].value_counts())

palette={
    'iPT ME 1': '#A7C8EB',
    'iPT ME 2': '#8CC47A',
    'iPT ME 3': '#C9AFD5',
    'iPT ME 4': '#F18C93',
    'iTAL ME 1': '#F0B836',
    'iTAL ME 2': '#007ABA',
    'iTAL ME 3':'#8F65A8',
    'iTAL ME 4': '#D975AC',
    'outlier':'#A8A8A7'
}##'#92D2DF',

color_dict = {
    'iPT ME 4': '#F18C93', # red
    'outlier': '#808080',   # black
    'iTAL ME 1':'#F0B836'
}

sc.pl.umap(adata3, color = 'for_plotting', palette=color_dict)
sc.pl.umap(adata3, color = 'for_plotting_iTAL', palette=color_dict)

sc.pl.scatter(
            adata3,
            x="CenterX_global_px",
            y="CenterY_global_px",
            color='for_plotting',
            size=3,
            title='whatever',
            frameon=False,
            color_map='plasma',
            legend_loc='on data',
            # Set the colormap to plasma
            # Adjusting the figure background color is not directly supported in sc.pl.scatter,
)



sc.pl.scatter(
            adata3,
            x="CenterX_global_px",
            y="CenterY_global_px",
            color='fov',
            size=3,
            title='whatever',
            frameon=False,
            color_map='plasma',
            legend_loc='on data',
            # Set the colormap to plasma
            # Adjusting the figure background color is not directly supported in sc.pl.scatter,
)

original_palette = {
    'PC': (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    'CNT': (0.6823529411764706, 0.7803921568627451, 0.9098039215686274),
    'DCT': (1.0, 0.4980392156862745, 0.054901960784313725),
    'DTL_ATL': (1.0, 0.7333333333333333, 0.47058823529411764),
    'EC_DVR': (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    'EC_Peritub': (0.596078431372549, 0.8745098039215686, 0.5411764705882353),
    'EC_glom': (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
    'IC A': (1.0, 0.596078431372549, 0.5882352941176471),
    'IC B': (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
    'Immune': (0.7725490196078432, 0.6901960784313725, 0.8352941176470589),
    'Podo': (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
    'Fibroblast': (0.7686274509803922, 0.611764705882353, 0.5803921568627451),
    'PEC': (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
    'PT': (0.9686274509803922, 0.7137254901960784, 0.8235294117647058),
    'MC1': (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
    'iPT': (0.7803921568627451, 0.7803921568627451, 0.7803921568627451),
    'iTAL': (0.8588235294117647, 0.8588235294117647, 0.5529411764705883),
    'TAL': (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
    'VSMC': (0.6196078431372549, 0.8549019607843137, 0.8980392156862745),
    'CD4+': (0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
    'Baso/Mast': (0.2235, 0.2314, 0.4745),
    'B': (0.7450980392156863, 0.7294117647058823, 0.8549019607843137),
    'CD8+ I': (0.6784, 0.2863, 0.2902),
    'Macro II': (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
    'Neutrophil': (0.9882352941176471, 0.803921568627451, 0.8980392156862745),
    'Macro I': (0.3, 0.3, 0.3),
    'NK': (0.7372549019607844, 0.5019607843137255, 0.7411764705882353),
    'cDC': (0.9176470588235294, 0.5019607843137255, 0.23529411764705882),
    'mDC': (0.34509803921568627, 0.33725490196078434, 0.6745098039215687),
    'pDC': (0.6509803921568628, 0.4627450980392157, 0.11372549019607843),
    'Macro III': (0.9294117647058824, 0.6941176470588235, 0.12549019607843137),
    'Macro IV': (0.41568627450980394, 0.23921568627450981, 0.6039215686274509),
    'cycling Lymphocytes': (0.8, 0.47058823529411764, 0.7372549019607844),
    'CD8+ II': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314),
    'Plasma': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
    'other':(0.5,0.5,0.5,0.5)
    #'cytotox B': (0.4, 0.7607843137254902, 0.6470588235294118)
}

##updated= iPT with color from DCT, Immune with color from glom endo but lighter, fibroblast with color of PC
palette = {
    'PC': (1, 1,1),
    'CNT': (1, 1,1),
    'DCT': (1, 1,1),
    'DTL_ATL': (1, 1,1),
    'EC_DVR': (1, 1,1),
    'EC_Peritub': (1, 1,1),
    'EC_glom': (1, 1,1),
    'IC A': (1, 1,1),
    'IC B': (1, 1,1),
    'Immune': (0.96, 0.57, 0.59),
    'Podo': (1, 1,1),
    'Fibroblast': (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    'PEC': (1, 1,1),
    'PT': (0.9686274509803922, 0.7137254901960784, 0.8235294117647058),
    'MC1': (1, 1,1),
    'iPT': (0.86, 0.78, 0.00),
    'iTAL': (0.8588235294117647, 0.8588235294117647, 0.5529411764705883),
    'TAL': (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
    'VSMC': (1, 1,1),
    'CD4+': (0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
    'Baso/Mast': (0.2235, 0.2314, 0.4745),
    'B': (0.7450980392156863, 0.7294117647058823, 0.8549019607843137),
    'CD8+ I': (0.6784, 0.2863, 0.2902),
    'Macro II': (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
    'Neutrophil': (0.9882352941176471, 0.803921568627451, 0.8980392156862745),
    'Macro I': (0.3, 0.3, 0.3),
    'NK': (0.7372549019607844, 0.5019607843137255, 0.7411764705882353),
    'cDC': (0.9176470588235294, 0.5019607843137255, 0.23529411764705882),
    'mDC': (0.34509803921568627, 0.33725490196078434, 0.6745098039215687),
    'pDC': (0.6509803921568628, 0.4627450980392157, 0.11372549019607843),
    'Macro III': (0.9294117647058824, 0.6941176470588235, 0.12549019607843137),
    'Macro IV': (0.41568627450980394, 0.23921568627450981, 0.6039215686274509),
    'cycling Lymphocytes': (0.8, 0.47058823529411764, 0.7372549019607844),
    'CD8+ II': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314),
    'Plasma': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
    'other':(0.5,0.5,0.5,0.5)
    #'cytotox B': (0.4, 0.7607843137254902, 0.6470588235294118)
}

adata3

adata= sc.read_h5ad('/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/raw_for_plotting/HK3035_raw.h5ad')

del adata.obsp['50_micron_connectivities']
del adata.obsp['50_micron_distances']

# Convert cell_id columns to sets
cell_ids_in_adata = set(adata.obs["cell_id"].to_list())
cell_ids_in_adata3 = set(adata3.obs["cell_id"].to_list())

# Find the intersection of cell IDs present in both objects
common_cell_ids = cell_ids_in_adata.intersection(cell_ids_in_adata3)

# Print the common cell IDs
#print(common_cell_ids)
mask = adata.obs["cell_id"].isin(common_cell_ids)
adata_subset = adata[mask]
adata_subset

del adata

columns_to_map=['for_plotting', 'for_plotting_iTAL']

for column in columns_to_map:

    #Create a dictionary from adata2 mapping cell_id to immune_cell_neighbor_calling
    immune_calling_map = adata3.obs.set_index('cell_id')[f'{column}'].to_dict()

    # Map immune_cell_neighbor_calling values to adata_subset based on cell_id
    adata_subset.obs[f'mapped_{column}'] = adata_subset.obs['cell_id'].map(immune_calling_map)

adata_subset

adata_subset.obsm['X_umap']=adata3.obsm['X_umap'].copy()
sc.pl.umap(adata_subset, color = "mapped_for_plotting", palette=color_dict)
sc.pl.umap(adata_subset, color = "mapped_for_plotting_iTAL", palette=color_dict)

adata_subset.obs['fov']=adata_subset.obs['fov'].astype(str)

import matplotlib.pyplot as plt
import squidpy as sq

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
#plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
#plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
#plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in adata.obs['fov'].unique():
        ax = sq.pl.spatial_segment(
            adata,
            color="mapped_for_plotting",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            vmax = 2,
            vmin = 0,
            colorbar=False,
            cmap='plasma',
            na_color=(0.3000000000001, 0.3000000000001, 0.3000000000001),

        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/iPT_subclustering/Images_3035/{fov}_iPTME4_HK3035.png'
        plt.savefig(filename, dpi=450, bbox_inches='tight', pad_inches=0)
        plt.show()

import matplotlib.pyplot as plt
import squidpy as sq

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
#plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
#plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
#plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in adata.obs['fov'].unique():
        ax = sq.pl.spatial_segment(
            adata,
            color="mapped_for_plotting_iTAL",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            vmax = 2,
            vmin = 0,
            colorbar=False,
            cmap='plasma',
            na_color=(0.3000000000001, 0.3000000000001, 0.3000000000001),

        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/iPT_subclustering/Images_3035/{fov}_iTALME1_HK3035.png'
        plt.savefig(filename, dpi=450, bbox_inches='tight', pad_inches=0)
        plt.show()

import matplotlib.pyplot as plt
import squidpy as sq

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
#plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
#plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
#plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]
fov_include_2 = ['3','4','5','11','12','13','19','20','21','27','28','29','35','36','37','43','44',
                 '45','50','51','52','56','57','58','60','61','62','67','68','69','75','76','77','83','84','85','91',
                 '92','93','99','100','101','107','108','109','115','116','117']

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in fov_include_2:
        ax = sq.pl.spatial_segment(
            adata,
            color="mapped_annotation_post_scanvi70_broad",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            vmax = 2,
            vmin = 0,
            colorbar=False,
            cmap='plasma',
            na_color=(0.3000000000001, 0.3000000000001, 0.3000000000001),

        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/iPT_subclustering/Images_2841/{fov}_HK2841.png'
        plt.savefig(filename, dpi=450, bbox_inches='tight', pad_inches=0)
        plt.show()

plt.style.use('dark_background')
adata_list = [adata_subset]
for adata in adata_list:
    adata_name = adata.obs["fov"].unique().astype(str)
    for fov in fov_include_2:
        ax = sq.pl.spatial_segment(
            adata,
            color="immune_organization_score",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            vmax = 1.3,
            vmin = 0,
            colorbar=True,
            cmap='plasma'
        )

        # Set the visibility of the spines to False
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/Figure_2/immune_black_1.3_organization_scalebar.png'

        # When saving, also specify the bbox_inches parameter to eliminate any extra white space
        plt.savefig(filename, dpi=900, bbox_inches='tight', pad_inches=0)
        plt.show()















adata_subset.obs["fov"]=adata_subset.obs["fov"].astype(str)

import matplotlib.pyplot as plt
import squidpy as sq
import gc
# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]
#fov_include_2 = ['1','2','3','11','12','13','22','23','24']

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in adata.obs['fov'].unique():
        ax = sq.pl.spatial_segment(
            adata,
            color="immune_organization_score",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            vmax = 2,
            vmin = 0,
            colorbar=False,
            cmap='plasma'
        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/Figure_2/immune_organization_{fov}_HK3039.png'
        plt.savefig(filename, dpi=900, transparent=True)
        plt.show()
        gc.collect()













import matplotlib.pyplot as plt
import squidpy as sq

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]
fov_include_2 = ['1','2','3','11','12','13','22','23','24']

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in adata.obs['fov'].unique():
        ax = sq.pl.spatial_segment(
            adata,
            color="mapped_annotation_post_scanvi70_broad",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            cmap='plasma',
            vmax = 2,
            vmin = 0,
            na_color=(0.2000000000001, 0.2000000000001, 0.2000000000001),
            scalebar_dx=0.12,
            scalebar_kwargs={"scale_loc": "bottom", "location": "upper right"},
        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'scale.png'
        plt.savefig(filename, dpi=900, transparent=True)
        plt.show()



adata_subset.obs["fov"]=adata_subset.obs["fov"].astype(str)

# Commented out IPython magic to ensure Python compatibility.
import matplotlib.pyplot as plt
import squidpy as sq

# %config InlineBackend.figure_format='retina'
plt.style.use('dark_background')

# Create your spatial segment plot
ax = sq.pl.spatial_segment(
    adata_subset,
    color="immune_coarse",
    library_key="fov",
    library_id=['1'],
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=0.12,
    scalebar_kwargs={"scale_loc": "bottom", "location": "upper right"},
    legend_fontsize='small',
    return_ax=True
)

# Adjust the current figure to make space for the legend below
plt.gcf().subplots_adjust(bottom=1.2)  # You may need to adjust this value

# Get the handles and labels from the plot
handles, labels = ax.get_legend_handles_labels()

# Create the legend under the plot
ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=4, frameon=False)  # Adjust ncol as needed
plt.savefig('TLO_legend_immune.png', dpi=600)
plt.show()

sc.pp.normalize_total(adata_subset, inplace=True)
sc.pp.log1p(adata_subset)

from scipy.sparse import issparse
markers = ["LTB"]
markers = [gene for gene in markers if gene in adata_subset.var_names]

# Subset the .X matrix for the markers and convert to a DataFrame
gene_expression_df = pd.DataFrame(
    adata_subset[:, markers].X.toarray() if issparse(adata_subset.X) else adata_subset[:, markers].X,
    index=adata_subset.obs_names,
    columns=markers
)
print(gene_expression_df.head(20))

from sklearn.preprocessing import MinMaxScaler

# Initialize the MinMaxScaler with the desired feature range
scaler = MinMaxScaler(feature_range=(0, 1))

# Select the 'PC' column and scale it
# Note: .values.reshape(-1, 1) converts it from 1D array to 2D array as expected by the scaler
pc_scaled = scaler.fit_transform(gene_expression_df[['LTB']].values.reshape(-1, 1))

# Replace the original 'PC' column with the scaled values
gene_expression_df['LTB_scaled'] = pc_scaled.flatten()
print(gene_expression_df.head(20))

adata_subset.obs["LTB_scaled"]=gene_expression_df['LTB_scaled']

adata_subset.obs["fov"]=adata_subset.obs["fov"].astype(str)

TNFSF13B_scaled

import matplotlib.pyplot as plt
import squidpy as sq

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]
fov_include_2 = ['1','2','3','11','12','13','22','23','24']

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in adata.obs['fov'].unique():
        ax = sq.pl.spatial_segment(
            adata,
            color="TNFSF13B_scaled",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            na_color=(0.2000000000001, 0.2000000000001, 0.2000000000001),
            cmap='plasma',
            colorbar=False,

        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/Immune_TLO/Penumbra_TLO/Pictures/TNFSF13B_{fov}_HK2695.png'
        plt.savefig(filename, dpi=900, transparent=True)
        plt.show()

import matplotlib.pyplot as plt
import squidpy as sq

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]
fov_include_2 = ['1','2','3','11','12','13','22','23','24']

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in adata.obs['fov'].unique():
        ax = sq.pl.spatial_segment(
            adata,
            color="TNFRSF13C_scaled",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            na_color=(0.2000000000001, 0.2000000000001, 0.2000000000001),
            cmap='plasma',
            colorbar=False,

        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/Immune_TLO/Penumbra_TLO/Pictures/TNFRSF13C_{fov}_HK2695.png'
        plt.savefig(filename, dpi=900, transparent=True)
        plt.show()

from scipy.sparse import issparse
markers = ["CCL21"]
markers = [gene for gene in markers if gene in adata_subset.var_names]

# Subset the .X matrix for the markers and convert to a DataFrame
gene_expression_df = pd.DataFrame(
    adata_subset[:, markers].X.toarray() if issparse(adata_subset.X) else adata_subset[:, markers].X,
    index=adata_subset.obs_names,
    columns=markers
)
print(gene_expression_df.head(20))

from sklearn.preprocessing import MinMaxScaler

# Initialize the MinMaxScaler with the desired feature range
scaler = MinMaxScaler(feature_range=(0, 1))

# Select the 'PC' column and scale it
# Note: .values.reshape(-1, 1) converts it from 1D array to 2D array as expected by the scaler
pc_scaled = scaler.fit_transform(gene_expression_df[['CCL21']].values.reshape(-1, 1))

# Replace the original 'PC' column with the scaled values
gene_expression_df['CCL21_scaled'] = pc_scaled.flatten()
print(gene_expression_df.head(20))

adata_subset.obs["CCL21_scaled"]=gene_expression_df['CCL21_scaled']

import matplotlib.pyplot as plt
import squidpy as sq

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures

adata_list = [adata_subset]
fov_include_2 = ['1','2','3','11','12','13','22','23','24']

for adata in adata_list:
    adata_name=adata.obs["fov"].unique().astype(str)
    for fov in adata.obs['fov'].unique():
        ax = sq.pl.spatial_segment(
            adata,
            color="CCL21_scaled",
            library_key="fov",
            library_id=fov,
            seg_cell_id="cell_ID",
            seg_outline=True,
            img=False,
            title='',
            axis_label='',
            return_ax=True,
            frameon=False,
            na_color=(0.2000000000001, 0.2000000000001, 0.2000000000001),
            cmap='plasma',
            colorbar=False,

        )

        # Remove the legend, if present
        if ax.get_legend():
            ax.get_legend().remove()

        # Corrected filename construction with string formatting
        filename = f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/Immune_TLO/Penumbra_TLO/Pictures/CCL21_{fov}_HK2695.png'
        plt.savefig(filename, dpi=900, transparent=True)
        plt.show()





import matplotlib.pyplot as plt
import squidpy as sq
import gc

# Apply dark background style
plt.style.use('dark_background')

# Manually adjust the background color to be transparent for subsequent plots
plt.rcParams['figure.facecolor'] = 'none'  # For the figure background
plt.rcParams['axes.facecolor'] = 'none'  # For the axes background
plt.rcParams['savefig.facecolor'] = 'none'  # For the saved figures


for library_id in adata_subset.obs["fov"].unique():
    print(library_id)
    ax = sq.pl.spatial_segment(
        adata_subset,
        color="CCL19_scaled",
        library_key="fov",
        library_id=[library_id],
        seg_cell_id="cell_ID",
        seg_outline=True,
        img=False,
        title='',
        axis_label='',
        return_ax=True,
        frameon=False,
        na_color=(0.2000000000001, 0.2000000000001, 0.2000000000001),
        cmap='plasma',
        colorbar=False
    )

    # Remove the legend, if present
    if ax.get_legend():
        ax.get_legend().remove()
    plt.savefig(f'/content/drive/MyDrive/Bernhard/Cosmx/After_Cleaning/Spatial_Neighborhoods/ACSM2B_scaled_{library_id}', dpi=600)
    plt.show()
    gc.collect()



adata_subset.obs["fov"]

