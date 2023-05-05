# my Python cheatsheet

## JuPyteR notebook

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


## os


