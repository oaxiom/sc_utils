# sc_utils
A series of utilities for scRNA-seq analysis

Also includes a memory safe(r) version of scran compute_sum_factors() that uses just-in-time
conversions to dense arrays and uses sparse arrays wherever possible.

To use:

```
adata = sc.read('./data.h5ad')

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added='quick_clusters', resolution=0.5, flavor='igraph', n_iterations=2)

compute_sum_factors(adata, clusters='quick_clusters',
                            parallelize=False, # True == broken?
                            algorithm='CVXPY', # Can use sparse matrices...
                            max_size=3000,
                            min_mean=0.1,
                            plotting=True,
                            lower_bound=0.4,
                            normalize_counts=False,
                            save_plots_dir='./scranpy/')
```
