# my Python cheatsheet

## Basics

### Lists

List comprehension:

```
[len(i) for i in mylist]
```

Search pattern in a list:

```
# search in list of strings
[i for i in mylist if 'ENSG' in i]

# here, search in AnnData.var index
[i for i in concatenatedData.var.index.tolist() if i.startswith('ENSG')]
[i for i in concatenatedData.var.index.tolist() if i.endswith('ENSG')]
```

Convert list of bytes to strings:

```
[i.decode("utf-8") for i in mylist]
```

### Sets

Create sets (unique, unordered):

```
a=set[1, 2, 3])
b=set[2, 1, 6]
```

```
# Intersection:
a.intersection(b)
a & b

# Union:
a.union(b)
a | b

# Difference:
a.difference(b)
a - b

# Symmetric Difference (elements from both sets that are not present on the other):
a.symmetric_difference(b)
a ^ b
```

## Jupyter notebook

Run Jupyter notebook from command line:

```
pip install runipy

runipy MyNotebook.ipynb
```

Install Python packages with pip from notebook block:

```
!pip install pandas
!pip install pandas --user
```

### os

Join paths (handles slashes properly)

```
os.path.join(mydirectory, "filename.txt")
```

### rpy2

Load:

```
%load_ext rpy2.ipython
```

Pass Python arguments (if run inside a notebook with Python kernel activated)

```
%%R -i myargument
```


## anndata

Read h5ad files into list and concatenate:

```
sampleData = []
for path in allh5ad:
    temp=sc.read_h5ad(str(path).strip())
    temp.obs_names_make_unique()
    sampleData.append(temp)
    
concatenatedData = anndata.concat(adatas=sampleData,
                                               keys = sampleIds, # should be unique
                                               join = "outer",
                                               axis=0,
                                               index_unique=None)
```

Write h5ad:

```
anndata.AnnData.write_h5ad(myadata, "adata.h5ad")
```

## scipy

Read a CSR (compressed sparse row) sparse matrix:

```
counts = sp.sparse.csr_matrix((data, indices, indptr), shape=(len(h5_features), len(h5_barcodes)))
```

Convert to CSC (compressed sparse column):

```
counts.tocsc()
```

Convert to dense format:

```
sp.sparse.csr_matrix.todense(counts)
```

## pandas

Order pandas data frame using another list:

```
df['A'] = pd.Categorical(df['A'], ordered=True, categories=list_order)
df = df.sort_values('A')
```

Add prefix to column names:

```
df=pd.DataFrame({'A': [1, 2, 3, 4], 'B': [3, 4, 5, 6]})
df.add_prefix('col_')
```

## scanpy

Multiple UMAPs + figure size control:

```
with rc_context({'figure.figsize': (8, 6)}):
    sc.pl.umap(adataMix2, 
           color=['Variable1', 'Variable2'],
          s=20, ncols=2, wspace=0.3)
```
Save figures:

```
sc.settings.figdir = os.path.join("dir_figures")


```

## h5py

Check content of hdf5 file:

```
out_h5="counts.hdf5"

with h5py.File(out_h5, "r") as check:
    h5_barcodes = list(check['barcodes'])
```

## matplotlib

```
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
```

### Figure control

```
plt.rcParams['figure.figsize'] = [12, 12]
```
### Barplots

```
mydf_counts=mydf.value_counts()
cluster = list(mydf_counts.keys())
cells = list(mydf_counts)
      
fig = plt.figure(figsize = (10, 5))
     
# creating the bar plot
plt.bar(cluster, cells, color ='maroon', width = 0.4)
     
plt.xlabel("Cluster name")
plt.ylabel("Number of cells")
plt.xticks(rotation=90)
plt.title("title")
plt.savefig("barplot.png)
plt.show()
```

### Venn diagrams

Example of a 3-way Venn

```
venn3(subsets=[set(list1), set(list2), set(list3)],
     set_labels=("list 1", "list 2", "list 3"))

plt.savefig('Venn3.png')
```

## Conversion of R / Seurat object into anndata / h5ad

```
library(sceasy)
obj <- readRDS("in_seurat.rds")
sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
                       outFile='out_anndata.h5ad')
```

Sometimes, using a "SingleCellExperiment" intermediate object can help with the conversion:

```
rds1 <- "in_seurat.rds"
sce1 <- "out_sce.rds"
h5ad1 <- "out_anndata.h5ad"

obj <- readRDS(rds1)

# Seurat to SCE
sceasy::convertFormat(obj, from="seurat", to="sce",
                        outFile=sce1,
                        assay="RNA",
                        main_layer="counts")
sce_in <- readRDS(sce1)

# SCE to AnnData
sceasy::convertFormat(sce_in, from="sce", to="anndata",
                       outFile=h5ad1)
```
