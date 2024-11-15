## Set up virtual environment
#python3.9 -m venv .myscenv
#source .myscenv/bin/activate # python version 3.9.18
## install packages
#python3 -m pip install --upgrade pip 
#python3 -m pip install -r pip_packages.txt

import os
import scanpy as sc
import pandas as pd
import seaborn as sns
import mudata as md
import anndata as ad
from datetime import datetime
import itertools
import re

# function to print out date and time
def timenow():
    return(datetime.now().strftime("%y-%d-%m %H:%M:%S"))

################ Set up ########################

## name of directory that will be created locally to output results and intermediate files
dir_name="sc_out/"

# name of local directory 
local_path="/home/"
os.makedirs(local_path, exist_ok=True)

# name of directory where to store resulting pre-processing intermediate files
local_out=os.path.join(local_path, dir_name)
os.makedirs(local_out, exist_ok=True)

# setting scanpy outpu figure directory
figures_out=os.path.join(local_out, "figures")
os.makedirs(figures_out, exist_ok=True)
sc.settings.figdir = figures_out

################## Sync and read in AnnData ##############################

## Input h5mu name
in_mu="myadata.h5mu"

# read in file
adata=md.read_h5mu(os.path.join(local_path, in_mu))["rna"]
adata.shape

# Make sure .X is raw
print(adata.X.min())
print(adata.X.max())

from scanpy._utils import check_nonnegative_integers
check_nonnegative_integers(adata.X)

################## Filter #####################

## Remove lowly expressed genes
sc.pp.filter_genes(adata, min_cells=3, inplace=True)
adata.shape

############# Normalization and log scaling #################

# Saving raw count data in layer "counts"
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)

# Logarithmize the data
sc.pp.log1p(adata)

# save log normalized data in layer "log_normalized"
adata.layers["log_normalized"]=adata.X.copy()

######### Calculate QC stats (will use total counts and %mito content) before selecting protein coding genes only #########

# mitochondrial genes
adata.var['mito'] = adata.var.index.str.startswith('MT-')
adata.var['mito'].value_counts()

# calculate QC metrics for mitochondrial and ribosomal genes
#sc.pp.calculate_qc_metrics(adata, qc_vars='mito', percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(adata, qc_vars='mito', percent_top=None, log1p=True, inplace=True)
    # 'n_genes_by_counts': number of genes with at least 1 count in a cell. Seurat's "nFeature_RNA" = nGene per cell
    # 'total_counts': sum of counts for a gene. Seurat's "nCount_RNA" = nUMI per gene
    # 'total_counts_mito': sum of counts for a cell which are mitochondrial
    # 'total_counts_ribo': --- ribosomal
    # 'pct_counts_mito': proportion of total counts for a cell which are mitochondrial.
    # 'pct_counts_ribo': --- ribosomal

########### Select protein coding genes ############

## Source (we could also fetch a gencode annotation file): https://osf.io/mhda7/ -> Human Protein-Coding genes in 2029 (downloaded Genes.xlsx and extract gene symbols as csv file)
protein_coding=pd.read_csv(os.path.join(local_path, "Genes_proteincoding_symbols_from_osf.csv"))["Gene_Symbol"].tolist()
len(protein_coding) 

# Extract genes from protein_coding that are present in adata.var
protein_coding_ls=list(set(adata.var.index) & set(protein_coding))
len(protein_coding_ls)

# subset to protein-coding genes only
adata=adata[:,protein_coding_ls]
adata.shape

############### Select HVG per batch and take union #################

# create empty lists
hvg_list3=[]
hvg_list4=[]
# create empty dictionaries
hvg_dic3={}
hvg_dic4={}

for i in adata.obs.Dataset.unique():
    print(i)
    # subset adata
    tmp=adata[adata.obs.Dataset==i]
    print(tmp.shape)
    # using top 2000 HVGs
    sc.pp.highly_variable_genes(tmp, batch_key="Dataset", n_top_genes=2000, subset=False)
    tmp_hvg=tmp.var.loc[tmp.var['highly_variable'],].index.tolist()
    print(len(tmp_hvg))
    hvg_list3.append(tmp_hvg)
    hvg_dic3[i]=tmp_hvg
    tmp.var.drop(columns=["highly_variable", "highly_variable_nbatches", "highly_variable_intersection"], inplace=True)
    # using mean, disp criteria
    sc.pp.highly_variable_genes(tmp, batch_key="Dataset", min_mean=0.0125, max_mean=3, min_disp=0.05, subset=False)
    tmp_hvg=tmp.var.loc[tmp.var['highly_variable'],].index.tolist()
    print(len(tmp_hvg))
    hvg_list4.append(tmp_hvg)
    hvg_dic4[i]=tmp_hvg
    
# make sure the list length is correct
len(hvg_list4)

# flatten lists
flat_hvg_list4 = list(itertools.chain(*hvg_list4))
unique_flat_hvg_list4=list(set(flat_hvg_list4)) # 11520

flat_hvg_list3 = list(itertools.chain(*hvg_list3))
unique_flat_hvg_list3=list(set(flat_hvg_list3)) # 11764

df_unique_flat_hvg_list4=pd.DataFrame(unique_flat_hvg_list4, columns=["HVG"])
df_unique_flat_hvg_list3=pd.DataFrame(unique_flat_hvg_list3, columns=["HVG"])

# intersection
len(set(unique_flat_hvg_list3).intersection(set(unique_flat_hvg_list4))) # 6158

df_unique_flat_hvg_list4.to_csv(os.path.join(local_out, pref_pre_regress+"_HVG_list_disp_mean_union.tsv"), sep="\t", index=False)
df_unique_flat_hvg_list3.to_csv(os.path.join(local_out, pref_pre_regress+"_HVG_list_top2000_union.tsv"), sep="\t", index=False)

#### store results in dataframe, keeping the number and name of datasets in which each gene is found HV
df_hvg_dic4=pd.DataFrame(columns=["HVG", "n_datasets", "Dataset"], index=range(nhvg))
df_hvg_dic4["HVG"]=unique_flat_hvg_list4

# loop through unique HVG
for i in unique_flat_hvg_list4:
    # set counter for # of datasets in which a gene is HV
    gcounter=0
    # set list for names of datasets in which a gene is HV
    lsdatasets=[]
    # loop through HVG dictionary
    for j in hvg_dic4:
        if i in hvg_dic4[j]:
            # if gene is HVG in that dataset, increase counter and add name to list
            gcounter=gcounter+1
            lsdatasets.append(j)
    #print(i + " is in " + str(gcounter) + " datasets")
    # fill up data frame
    df_hvg_dic4.loc[df_hvg_dic4.HVG==i,["n_datasets", "Dataset"]]=[gcounter, ','.join(lsdatasets)]

print(df_hvg_dic4.head())
print(df_hvg_dic4.tail())

## Save to file
df_hvg_dic4.sort_values("n_datasets", ascending=False).to_csv(os.path.join(local_out, pref_pre_regress+"_DF_HVG.tsv"), sep="\t", index=False)

## Select genes that are highly variable in AT LEAST 2 datasets
ls_hvg_selection=df_hvg_dic4.loc[df_hvg_dic4.n_datasets >= 2, "HVG"].tolist()
nhvg=len(ls_hvg_selection)

# prefix to read/save files
pref_pre_regress="proteincoding_"+str(nhvg)+"HVG"

## selection of these nhvg genes in flat_hvg_list4
ls_hvg_selection_flat=[i for i in flat_hvg_list4 if i in ls_hvg_selection]

## Prep new columns for "adata.var"
# data frame from full gene list +  count in how many datasets each gene is HVG
var_new=pd.DataFrame({"feature_name":ls_hvg_selection_flat}).groupby("feature_name").size()
# create a data frame
var_new2=pd.DataFrame({"highly_variable_nbatches":var_new, "highly_variable": True})
var_new2.shape

## Modify adata.var and save file
adata.var=adata.var.merge(var_new2, how="left", left_index=True, right_index=True)
adata.var.loc[adata.var.highly_variable.isnull(), ["highly_variable_nbatches", "highly_variable"]]=[0, False]
adata.var.head()

# check if the process was done correctly
adata.var.loc[adata.var.highly_variable==True,].highly_variable_nbatches.sort_values(ascending=False)

# convert to boolean to be able to use the column for subsetting easily
adata.var.highly_variable=adata.var.highly_variable.astype('bool')

md.MuData({"rna": adata}).write(os.path.join(local_out, pref_pre_regress+"_norm_hvg.h5mu"), compression="gzip")

################ Regressing out unwanted variability ####################

## Cell cycle genes
g1s_genes = [x.strip() for x in open(os.path.join(local_path, "G1S_genes.txt"))]
g1s_genes = [x for x in g1s_genes if x in adata.var_names]

g2m_genes = [x.strip() for x in open(os.path.join(local_path, "G2M_genes.txt"))]
g2m_genes = [x for x in g2m_genes if x in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=g1s_genes, g2m_genes=g2m_genes, copy=False)

# check variables
adata.obs[['total_counts', 'pct_counts_mito', "S_score", "G2M_score"]].head(3)

# save up to here as raw (before subsetting to HVG only). 
# The .X from raw contains log normalized counts for all genes.
adata.raw = adata.copy()

# save scores for all genes
md.MuData({"rna": adata}).write(os.path.join(local_out, pref_pre_regress+"_scoresregression.h5mu"), compression="gzip")

### Subset to HVG only for regression ###
adata = adata[:, adata.var.highly_variable]
adata.shape

##################################################################
############# NO regression: scale, pca, nn, umap ################
##################################################################

adata_backup=adata.copy()

# scale
sc.pp.scale(adata, max_value=10)
md.MuData({"rna": adata}).write(os.path.join(local_out, pref_pre_regress+"_scaled.h5mu"), compression="gzip")
## PCA, NN, UMAP
sc.tl.pca(adata, svd_solver="arpack")
md.MuData({"rna": adata}).write(os.path.join(local_out, pref_pre_regress+"_pca.h5mu"), compression="gzip")

for i in ["Dataset", "CellType"]:
    sc.pl.pca(adata, color=i, save="_"+pref_pre_regress+"_"+i+".png")

# NN + UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

## Plots
for i in ["Dataset", "CellType"]:
    sc.pl.umap(adata, color=i, save="_"+pref_pre_regress+"_"+i+".png")

################################################################################
##################### Regress out, scale, pca, nn, umap ########################
################################################################################

# restore HVG-selected AnnData without PCA, UMAP etc.
adata=adata_backup.copy()

## new prefix post regression
pref="Regress_proteincoding_nointeg_"+str(nhvg)+"HVG"

# regress out 3 variables: total counts, mitochondrial content and cell cycle genes
print(timenow())
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mito", "S_score", "G2M_score"])
print(timenow())

## scale
sc.pp.scale(adata, max_value=10)

## PCA
# what is given to harmony are the PC coordinates, calculated based on the regressed normalized data.
sc.tl.pca(adata, svd_solver="arpack")
md.MuData({"rna": adata}).write(os.path.join(local_out, pref+"_pca.h5mu"), compression="gzip")

for i in ["Dataset", "CellType"]:
    sc.pl.pca(adata, color=i, save="_"+pref+"_"+i+".png")

# NN + UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)
md.MuData({"rna": adata}).write(os.path.join(local_out, pref+"_umap.h5mu"), compression="gzip")

## Plots
for i in ["Dataset", "CellType"]:
    sc.pl.umap(adata, color=i, save="_"+pref+"_"+i+".png")

############# pip freeze ################

# Save pip freeze / package versions
cmd="pip freeze > " + os.path.join(local_out, "pip_freeze.txt")
os.system(cmd)
