
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

