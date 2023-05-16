# my Python cheatsheet

## Basics

### Lists

List comprehension:

```
[len(i) for i in mylist]
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


## scanpy

Multiple UMAPs + figure size control:

```
with rc_context({'figure.figsize': (8, 6)}):
    sc.pl.umap(adataMix2, 
           color=['Variable1', 'Variable2'],
          s=20, ncols=2, wspace=0.3)
```

## h5py

Check content of hdf5 file:

```
out_h5="counts.hdf5"

with h5py.File(out_h5, "r") as check:
    h5_barcodes = list(check['barcodes'])
```

