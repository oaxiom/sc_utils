
import pandas as pd
from anndata import AnnData

def sparsify(filename):
    '''
    **Purpose**
        Convert a dense array in filename into a sparse array and return a
        AnnData object

    '''
    print('Sparsifying {0}'.format(filename))
    data = pd.read_csv(filename, index_col=0, header=0)
    genes = data.columns
    cells = data.index
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')

    '''
    oh = open('gene_names.{0}.tsv'.format(os.path.split(filename)[1]), 'w')
    for g in genes:
        oh.write('%s\n' % g)
    oh.close()
    '''

    print('Loaded {0}'.format(filename))
    ad = AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data
    return ad

def export_dense(adata, gene_name_filename, group_filename, dense_filename):
    '''
    **Purpose**
        Export a dense representation for bug testing and things like scran and
        seurat
    '''
    oh = open(gene_name_filename, 'w')
    for g in adata.var_names:
        oh.write('%s\n' % g)
    oh.close()

    # Save the groups:
    print('Save groups')
    input_groups = adata_pp.obs['groups']
    input_groups.to_csv(group_filename, sep='\t')
    print(input_groups)
    del adata_pp

    print('Save dense matrix')
    # save data_mat and go into R:
    #data_mat = adata.X.T.tocsr() # This one is dense already
    data_mat = adata.X # This one is dense already
    # Save out a dense array of the sparse array:
    oh = open(dense_filename, 'w')
    for i in range(data_mat.shape[0]):

        if i % 100 == 0:
            print('%s/%s' % (i, data_mat.shape[0]))

        oh.write('\t'.join([str(i) for i in data_mat.getrow(i).toarray()[0,:]]))
    oh.close()
