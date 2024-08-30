# -*- coding: utf-8 -*-
"""Immune Cosmx/SN/SC Atlas SCVI SCANVI.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/13o0D2AXReMKv8Ds38hf_XUUoVHLbkBuV
"""

!pip install matplotlib==3.5.3

!pip install torch==1.13.1
!pip install torchaudio==0.13.1
!pip install torchmetrics==0.11.4
!pip install torchsummary==1.5.1
!pip install torchtext==0.14.1
!pip install torchvision==0.14.1

# Commented out IPython magic to ensure Python compatibility.
#restart session

!pip install scanpy==1.8.2
!pip install rich==13.2


!pip install --quiet scvi-colab==0.11.0
from scvi_colab import install
install()


import rich

import sys
IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    !pip install --quiet scrublet
    import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import anndata
import scvi
import scanpy as sc

sc.set_figure_params(figsize=(4, 4))
scvi.settings.seed = 94705

# %config InlineBackend.print_figure_kwargs={'facecolor' : "w"}
# %config InlineBackend.figure_format='retina'

!pip install anndata2ri
import anndata2ri

!pip install session_info

seed=10
#sc.logging.print_versions()

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SN_SC_Cosmx_Immune_preSCVI.h5ad')
adata

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="tech",
    layer="counts", categorical_covariate_keys=["proj", "orig_ident"],
    continuous_covariate_keys=["nCount_RNA", "percent_mt"])
model = scvi.model.SCVI(adata)
model
vae = scvi.model.SCVI(adata, n_layers=5, n_latent=30, gene_likelihood="nb", dropout_rate=0.1)
vae.train(max_epochs = 600, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 25)
model = vae
model.save("/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SCVI_Models/Model_SN_SC_Cosmx_1", overwrite = True)

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4)

sc.pp.neighbors(adata, n_pcs=30, use_rep="X_scVI", random_state=seed)
sc.tl.umap(adata, min_dist=0.3, random_state=seed)

#sc.tl.leiden(adata, key_added="leiden_scVI_0_5", resolution=0.5, random_state=seed)
sc.tl.leiden(adata, key_added="leiden_scVI_0_3", resolution=0.3, random_state=seed)
sc.tl.leiden(adata, key_added="leiden_scVI_0_1", resolution=0.1, random_state=seed)

adata.write('/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SN_SC_Cosmx_Immune_SCVI_V1.h5ad')

#adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SN_SC_Cosmx_Immune_SCVI_V1_SCANVI_annotated_V1.h5ad')

adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SN_SC_Cosmx_Immune_SCVI_V1_SCANVI_annotated_V5.h5ad')

#adata=sc.read_h5ad('/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SN_SC_Cosmx_Immune_SCVI_V1_SCANVI_V2.h5ad')

print(adata.obs["annotations_after_scanvi"].value_counts())

model = scvi.model.SCVI.load("/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SCVI_Models/Model_SN_SC_Cosmx_1", adata=adata)

seed = 10
scvi.settings.seed = 10

model_path = "/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SCVI_Models/Model_SN_SC_Cosmx_annotated_SCANVI_5"

scanvi_model = scvi.model.SCANVI.from_scvi_model(model, labels_key="annotations_after_scanvi", unlabeled_category = "Unknown")

scanvi_model.train(max_epochs = 600, plan_kwargs={"lr":0.001}, early_stopping = True, early_stopping_patience = 15)

scanvi_model.save(model_path, overwrite=True)

adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

sc.tl.leiden(adata, key_added="leiden_scanvi_0_5", resolution=0.5, random_state=seed)
sc.tl.leiden(adata, key_added="leiden_scanvi_0_3", resolution=0.3, random_state=seed)
sc.tl.leiden(adata, key_added="leiden_scanvi_0_1", resolution=0.1, random_state=seed)

adata.write('/content/drive/MyDrive/Bernhard/Immune_Cell_Atlas/SN_SC_Cosmx_Immune_SCVI_V1_SCANVI_annotated_V6.h5ad')