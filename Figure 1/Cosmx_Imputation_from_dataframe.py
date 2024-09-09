#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
import leidenalg as la
from pathlib import Path
from anndata import AnnData


# In[2]:


print('hallo')


# In[3]:


adata = sc.read_h5ad('/home/bcd/cosmx/All_Runs/Cosmx_SN_cleaned_inner_join.h5ad')


# In[3]:


#adata


# In[4]:


SN = sc.read_h5ad('/home/bcd/Atlas_sn/Atlas6_Human.h5ad')


# In[5]:


# Get the list of cell IDs from both objects
cell_ids_in_adata = adata.obs.index
cell_ids_in_sn = SN.obs.index

# Find the intersection of cell IDs present in both objects
common_cell_ids = cell_ids_in_adata.intersection(cell_ids_in_sn)

# Subset the single_cell_object to include only the common cells
single_nuc = SN[common_cell_ids].copy()

# Now single_cell_object_subset contains only the cells that are also present in adata
single_nuc


# In[6]:


single_nuc.obs["tech"]= "SN_SEQ"


# In[7]:


single_nuc.obs["sample"]=single_nuc.obs["orig.ident"]


# In[8]:


del single_nuc.obsm["X_pca"]
del single_nuc.obsm["X_umap"]
single_nuc


# In[9]:


adata_cosmx= adata[adata.obs["tech"]=="Cosmx"]
del adata_cosmx.obsm["X_scANVI"]
del adata_cosmx.obsm["X_scVI"]
del adata_cosmx.obsm["X_umap"]
del adata_cosmx.obsm["_scvi_extra_categorical_covs"]
del adata_cosmx.obsm["_scvi_extra_continuous_covs"]
del adata_cosmx.layers["scvi_normalized"]
adata_cosmx


# In[10]:


adata_combined_var = ad.concat([adata_cosmx,single_nuc], join = "outer")


# In[11]:


keep_str = ["sample",'tech',"cell_ID"]
keep_numer = []
keep = keep_str #+ keep_numer
adata_combined_var.obs = adata_combined_var.obs[keep]

# ensure proper data types
for i in keep_str:
    adata_combined_var.obs[i]= adata_combined_var.obs[i].astype('str')

#for i in keep_numer:
 #   adata_combined_var.obs[i]= adata_combined_var.obs[i].astype('float')


# In[12]:


#adata_combined_var


# In[12]:


adata_combined_var = adata_combined_var[adata_combined_var.obs_names.sort_values(), :]
adata = adata[adata.obs_names.sort_values(), :]


# In[13]:


print(adata.obs["tech"])
print(adata_combined_var.obs["tech"])


# In[15]:


adata_CosMx = adata_combined_var[adata_combined_var.obs["tech"]=="Cosmx"]


# In[51]:


#adata_CosMx.write('/home/bcd/cosmx/All_Runs/Imputation/Cosmx_sn_imputed_all_genes_sparse_batches.h5ad')


# In[4]:


#adata_CosMx= sc.read_h5ad('/home/bcd/cosmx/All_Runs/Imputation/Cosmx_sn_imputed_all_genes_sparse_batches.h5ad')


# In[16]:


#from sklearn.neighbors import NearestNeighbors
#Cosmx_cells_mask = (adata.obs['tech'] == 'Cosmx')
#snRNA_cells_mask = (adata.obs['tech'] != 'Cosmx')

#CosMx_index = np.where(Cosmx_cells_mask)[0]
#snRNA_index = np.where(snRNA_cells_mask)[0]

#nn = NearestNeighbors(n_neighbors=15, metric='euclidean')
#nn.fit(adata.obsm["X_scANVI"][snRNA_cells_mask])

#distances_all_to_non_Cosmx, indices_all_to_non_Cosmx = nn.kneighbors(adata.obsm["X_scANVI"][Cosmx_cells_mask])


# In[14]:


distances_all_to_non_Cosmx = np.load('/home/bcd/cosmx/All_Runs/Imputation/distances_all_to_non_Cosmx.npy')
indices_all_to_non_Cosmx = np.load('/home/bcd/cosmx/All_Runs/Imputation/indices_all_to_non_Cosmx.npy')


# In[15]:


Cosmx_cells_mask = (adata.obs['tech'] == 'Cosmx')
snRNA_cells_mask = (adata.obs['tech'] != 'Cosmx')

CosMx_index = np.where(Cosmx_cells_mask)[0]
snRNA_index = np.where(snRNA_cells_mask)[0]


# In[16]:


adata_combined_var


# In[20]:


adata_CosMx = adata_combined_var[adata_combined_var.obs["tech"]=="Cosmx"]
adata_CosMx


# In[21]:


adata_CosMx.obs["cell_ID"]


# In[ ]:





# In[ ]:






# In[ ]:


import pandas as pd
from scipy.sparse import csr_matrix, vstack
import gc

# Specify the path to your CSV file
csv_file_path = '/home/bcd/cosmx/All_Runs/Imputation/Cosmx_all_genes.csv'

# Define the chunk size
chunk_size = 400000  # Number of rows per chunk

# Initialize a variable to hold the final CSR matrix
final_csr_matrix = None

# Create an iterator object for reading in chunks with the first column as the index
chunk_iter = pd.read_csv(csv_file_path, chunksize=chunk_size, index_col=0)

# Process the file in chunks
for chunk in chunk_iter:
    # Convert the current chunk to a CSR matrix
    current_csr_chunk = csr_matrix(chunk.values, dtype='float32')
    
    # If final_csr_matrix is None, this is the first chunk
    if final_csr_matrix is None:
        final_csr_matrix = current_csr_chunk
    else:
        # Stack the current chunk vertically with the previous chunks
        final_csr_matrix = vstack([final_csr_matrix, current_csr_chunk])
    print(chunk)
    gc.collect()





# In[ ]:


print(final_csr_matrix.shape)


# In[ ]:


adata_CosMx.layers["scanvi_imputed"] = final_csr_matrix
adata_CosMx.X = adata_CosMx.layers["scanvi_imputed"]


# In[ ]:


adata_CosMx.write('/home/bcd/cosmx/All_Runs/Imputation/Cosmx_sn_imputed_all_genes_from_dataframe.h5ad')
