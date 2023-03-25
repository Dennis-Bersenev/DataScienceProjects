#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

adata.var_names_make_unique()
adata


# In[2]:


adata.layers["counts"] = adata.X.copy()


# In[3]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# In[4]:


adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[5]:


adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]


# In[6]:


sc.pp.normalize_total(adata, target_sum=1e4)


# In[7]:


sc.pp.log1p(adata)


# In[8]:


sc.pp.highly_variable_genes(adata, min_mean=0.0005, max_mean=10, min_disp=0.1) 
# sc.pp.highly_variable_genes(adata, n_top_genes=4000) 


# In[9]:


adata.raw = adata
adata


# In[10]:


adata = adata[:, adata.var.highly_variable]


# In[11]:


sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])


# In[12]:


sc.pp.scale(adata, max_value=10)


# In[13]:


#PCA 
sc.tl.pca(adata, svd_solver='arpack')


# In[14]:


#Clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)


# In[15]:


sc.tl.umap(adata)


# In[16]:


sc.tl.leiden(adata, resolution=0.75, key_added="leiden_2")


# In[17]:


sc.tl.rank_genes_groups(adata, 'leiden_2', method='logreg')


# |Cluster ID	|Markers         |Cell Type       |
# |-----------|----------------|----------------|
# |0	        |IL7R,S100A4     |Memory CD4+     |
# |1	        |CD14,LYZ	     |CD14+ Mono      |
# |2	        |IL7R,CCR7	     |Naive CD4+ T    |
# |3	        |MS4A1	         |B               |
# |4	        |CD8A	         |CD8+ T          |
# |5	        |FCGR3A,MS4A7	 |FCGR3A+         |
# |6	        |GNLY,NKG7	     |NK              |
# |7	        |FCER1A,CST3     |Dendritic C     |
# |8	        |PPBP	         |Megakaryocytes  |
# 
# (cluster ID column is from Seurat, so in your case, match Markers with cell types, by getting the markers for each cluster!
#  The point of doing this is as a clustering sanity check)

# In[18]:


new_cluster_names = [
    'Memory CD4+', 'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'FCGR3A Monocytes', 'NK',
    'Dendritic', 'Megakaryocytes']
adata.rename_categories('leiden_2', new_cluster_names)


# In[19]:


pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)


# In[20]:


sc.pl.umap(adata, color='leiden_2', legend_loc='on data', title='', frameon=False, save='.pdf')


# In[21]:


def get_column(adata_instance, X, goi):
    goi_id = adata_instance.var.gene_ids[goi]
    goi_index = (adata_instance.var.gene_ids==goi_id).argmax()
    return X[:, goi_index] 

def highly_coexpressed(goi, cell_type=None):
    if cell_type == None:
        X = adata.X
    else:
        #run on subset
        X = adata.X[adata.obs.leiden_2 == cell_type]

    corr_X = np.corrcoef(np.transpose(X)) 
    goi_corr_col = get_column(adata, corr_X, goi)
    indices = np.argsort(goi_corr_col)
    n = indices.shape[0]
    indices_top_20 = np.flip(indices[n-21: n]) 
    temp = adata.var.gene_ids[indices_top_20].keys()
    d = {'Genes': adata.var.gene_ids[indices_top_20].keys(), 'Pearson Coefficients': goi_corr_col[indices_top_20]}
    df =  pd.DataFrame(data=d)
    return df

def pairwise_expression(goi_1, goi_2, cell_type=None):
    if cell_type == None:
        X = adata.X
    else:
        #run on subset
        X = adata.X[adata.obs.leiden_2 == cell_type]
    
    goi_1_col = get_column(adata, X, goi_1)
    goi_2_col = get_column(adata, X, goi_2)
    
    pw_corr = np.corrcoef(goi_1_col, goi_2_col)
    return pw_corr


# In[22]:


res = highly_coexpressed('LYZ')
res


# In[23]:


res2 = highly_coexpressed('LYZ', 'CD14 Monocytes')
res2


# In[24]:


pairwise_expression('LYZ', 'FTL') #It makes sense, it is the correlation between two vectors, outta be a number


# In[25]:


#Change adata to only include the monocyte rows!
pairwise_expression('LYZ', 'FTL', 'CD14 Monocytes')


# In[26]:


# adata.X = adata.layers["counts"].copy()
# sc.pl.scatter(adata, x='LYZ', y='FTL', use_raw=True) #Should recover OG plot, if using raw will plot log normalised values (which is what should be submitted!)


# In[27]:


# bdata = adata[adata.obs.leiden_2 == "CD14 Monocytes"] #Would be helpful to rename leiden_2 col to cell_type!
# sc.pl.scatter(bdata, x='LYZ', y='FTL', use_raw=True)

