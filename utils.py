
import sys, os, time, numpy
import gzip as gzipfile
import scipy as sp
import pandas as pd
from anndata import AnnData
import scanpy as sc
import matplotlib.cm as cm
import numpy as np
from matplotlib import colors
from .progress import progressbar

colors1 = cm.Greys_r(np.linspace(0.8,0.9,20))
colors2 = cm.Reds(np.linspace(0.0, 1, 100))
colorsComb = np.vstack([colors1, colors2])
cmap_grey_red = colors.LinearSegmentedColormap.from_list('cmap_grey_red', colorsComb)

def fastqPE(filename1, filename2, gzip=True):
    """
    generator object to parse fastQ PE files

    @HWI-M00955:51:000000000-A8WTD:1:1101:13770:1659 1:N:0:NNTNNNAGNNCNCTAT
    NGGTAAATGCGGGAGCTCCGCGCGCANNTGCGGCNNNGCATTGCCCATAATNNNNNNNCTACCGACGCTGACTNNNNNCTGTCTCTTATACACATNNNNGAGCCCACGNNNNCNNNCTAGNNNNNNNNNNNNNNNTTCTGCTTGTAAACA
    +
    #,,5,</<-<+++5+568A+6+5+++##5+5++5###+5+55-55A-A--5#######55+5<)+4)43++14#####*1*1*2011*0*1*1*1####***111(/'####/###-(((###############/-(/((./(((((((

    """
    if gzip:
        oh1 = gzipfile.open(filename1, "rt")
        oh2 = gzipfile.open(filename2, "rt")
    else:
        oh1 = open(filename1, "rt")
        oh2 = open(filename2, "rt")

    name1 = "dummy"
    while name1 != "":
        name1 = oh1.readline().strip()
        seq1 = oh1.readline().strip()
        strand1 = oh1.readline().strip()
        qual1 = oh1.readline().strip()

        name2 = oh2.readline().strip()
        seq2 = oh2.readline().strip()
        strand2 = oh2.readline().strip()
        qual2 = oh2.readline().strip()

        res = ({"name": name1, "strand": strand1, "seq": seq1, "qual": qual1},
            {"name": name2, "strand": strand2, "seq": seq2, "qual": qual2})
        yield res
    return

def __load_matrix(filename, csv=True):
    '''
    **Purpose**
        Convert a dense array in filename into a sparse array and return a
        AnnData object

    This is slower than just loading with Pandas!

    '''
    data = []
    barcodes = []
    if csv:
        sep = ','
    else:
        sep = '\t'

    oh = gzipfile.open(filename, 'rt')
    for i, line in enumerate(oh):
        if (i+1) % 10000 == 0:
            print(f'{i+1:,}')

        if i == 0:
            genes = line.split(sep)[1:]
            continue

        t = line.split(sep)
        cts = [float(c) for c in t[1:]]
        barcodes.append(t[0])
        data.append(cts)

    data = pd.DataFrame(data, index=barcodes, columns=genes) # i.e. transposed

    return data

def __load_genes_barcodes(filename, csv=True):
    '''
    **Purpose**
        Load in genes and barcodes

    '''
    barcodes = []
    if csv: sep = ','
    else: sep = '\t'

    oh = gzipfile.open(filename, 'rt')
    for i, line in enumerate(oh):
        if i == 0:
            genes = line.strip().split(sep)[1:]
            continue

        t = line.split(sep)
        barcodes.append(t[0])

    return barcodes, genes

def _drop_fusions_mir(data,
    gene_names,
    gene_ensg,
    ensg_to_symbol,
    drop_fusions:bool = True,
    drop_mir:bool = True,
    drop_tes:bool = True,
    drop_ribosomes:bool = False,
    ):

    # ensg_to_symbol must be valid to get here;
    todrop = []
    record_of_drops = []
    # drop genes with - in the form, but not -AS
    for n, e in zip(gene_names, gene_ensg):
        if '?' in e:
            todrop.append(e)
            record_of_drops.append(n)

        if drop_fusions and '-' in n:
            if n.startswith('MT-'): continue # Used for QC
            if n.startswith('ERV'): continue # Don't drop this ERVs, it's not a fusion!
            if 'ENS' not in e: continue
            if '-AS' in n: continue
            if '-int' in n: continue
            if 'Nkx' in n: continue # mouse genes
            if 'Krtap' in n: continue
            if ':' in n: continue # Don't drop TEs!
            todrop.append(e)
            record_of_drops.append(n)

        if drop_mir and n[0:3] == 'MIR' or n[0:3] == 'mir':
            if 'ENS' not in e: continue
            if 'hg' in e.lower(): continue # host genes;
            if ':' in n: continue # Don't drop TEs!
            todrop.append(e)
            record_of_drops.append(n)

        if drop_tes and ':' in n: # This only works with te_count;
            todrop.append(e)
            record_of_drops.append(n)

        if drop_ribosomes and n.startswith('RP'): # This only works with te_count;
            if n.startswith('RPL') or n.startswith('RPS'):
                todrop.append(e)
                record_of_drops.append(n)

        if drop_ribosomes and n.startswith('MRP'):
            if n.startswith('MRPL') or n.startswith('MRPS'):
                todrop.append(e)
                record_of_drops.append(n)

        if drop_ribosomes and n.startswith('Rp'): # This only works with te_count;
            if n.startswith('Rpl') or n.startswith('Rps'):
                todrop.append(e)
                record_of_drops.append(n)

    data.drop(todrop, axis=1, inplace=True)
    gene_names = []
    gene_ensg = data.columns # remap; # rebuild the gene names inde to avoid probelms with duplicate name/ensg? drops
    for ensg in gene_ensg:
        if ensg not in ensg_to_symbol:
            gene_names.append(ensg)
            #print(f'Warning: {ensg} not found')
        else:
            gene_names.append(ensg_to_symbol[ensg])
    return todrop, gene_names, gene_ensg

def _load_velocyte_mtx(path, load_ambiguous=False):

    ambiguous = None
    if load_ambiguous:
        print('Loading Spliced, Unspliced and Ambiguous matrices')
    else:
        print('Loading Spliced and Unspliced matrices')

    try:
        mtxU = np.loadtxt(os.path.join(path, 'Velocyto/raw/unspliced.mtx'), skiprows=3, delimiter=' ')
        shapeU = np.loadtxt(os.path.join(path, 'Velocyto/raw/unspliced.mtx'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)
    except FileNotFoundError:
        mtxU = np.loadtxt(os.path.join(path, 'Velocyto/raw/unspliced.mtx.gz'), skiprows=3, delimiter=' ')
        shapeU = np.loadtxt(os.path.join(path, 'Velocyto/raw/unspliced.mtx.gz'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)

    try:
        mtxS = np.loadtxt(os.path.join(path, 'Velocyto/raw/spliced.mtx'), skiprows=3, delimiter=' ')
        shapeS = np.loadtxt(os.path.join(path, 'Velocyto/raw/spliced.mtx'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)
    except FileNotFoundError:
        mtxS = np.loadtxt(os.path.join(path, 'Velocyto/raw/spliced.mtx.gz'), skiprows=3, delimiter=' ')
        shapeS = np.loadtxt(os.path.join(path, 'Velocyto/raw/spliced.mtx.gz'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects

    spliced = sp.sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
    unspliced = sp.sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()

    if load_ambiguous:
        try:
            mtxA = np.loadtxt(os.path.join(path, 'Velocyto/raw/ambiguous.mtx'), skiprows=3, delimiter=' ')
            shapeA = np.loadtxt(os.path.join(path, 'Velocyto/raw/ambiguous.mtx'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)
        except FileNotFoundError:
            mtxA = np.loadtxt(os.path.join(path, 'Velocyto/raw/ambiguous.mtx.gz'), skiprows=3, delimiter=' ')
            shapeA = np.loadtxt(os.path.join(path, 'Velocyto/raw/ambiguous.mtx.gz'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)

        ambiguous = sp.sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()

    try:
        genes = pd.read_csv(os.path.join(path, 'Velocyto/raw/features.tsv'), sep='\t',
            names=('gene_ids', 'feature_types'), index_col=1)
    except FileNotFoundError:
        genes = pd.read_csv(os.path.join(path, 'Velocyto/raw/features.tsv.gz'), sep='\t',
            names=('gene_ids', 'feature_types'), index_col=1)

    try:
        barcodes = pd.read_csv(os.path.join(path, 'Velocyto/raw/barcodes.tsv'), header=None, index_col=0)
    except:
        barcodes = pd.read_csv(os.path.join(path, 'Velocyto/raw/barcodes.tsv.gz'), header=None, index_col=0)

    barcodes.index.name = None

    return spliced, unspliced, ambiguous, genes, barcodes

def sparsify(filename=None, pandas_data_frame=None,
    obs_add=None,
    csv:bool = False,
    drop_fusions:bool = False,
    drop_mir:bool = False,
    drop_tes:bool = False,
    drop_ribosomes:bool = False,
    ensg_to_symbol = None,
    load_ambiguous:bool = False,
    velocyte_data = None):
    '''
    **Purpose**
        Convert a dense array in filename into a sparse array and return a
        AnnData object

    '''
    assert filename or (pandas_data_frame is not None), 'You must specify one of filename or pandas_data_frame'
    assert not (filename and (pandas_data_frame is not None)), 'You must specify only one of filename or pandas_data_frame'
    assert obs_add, 'obs_add needs to be provided'

    if drop_mir: assert ensg_to_symbol, 'Asked to drop_mir, but no ensg_to_symbol specified'
    if drop_fusions: assert ensg_to_symbol, 'Asked to drop_mir, but no ensg_to_symbol specified'

    print('Started {}'.format(filename))
    s = time.time()
    if filename:
        if csv:
            pdata = pd.read_csv(filename, index_col=0, header=0,
                sep=',',
                encoding='utf-8', compression='gzip', dtype={'x': int},
                low_memory=False,
                engine='c')

            # Numpy:
            barcodes, genes = __load_genes_barcodes(filename, csv=csv)
            npdata = np.loadtxt(filename, delimiter=',', skiprows=1, usecols=range(1, len(genes)+1))

            data = pd.DataFrame(npdata, index=barcodes, columns=genes)
        else:
            pdata = pd.read_csv(filename, index_col=0, header=0, sep='\t',
                encoding='utf-8', compression='gzip', dtype={'x': int},
                low_memory=False,
                engine='c')

            # Numpy:
            barcodes, genes = __load_genes_barcodes(filename, csv=csv)
            npdata = np.loadtxt(filename, delimiter='\t', skiprows=1, usecols=range(1, len(genes)+1))

            data = pd.DataFrame(npdata, index=barcodes, columns=genes)
    else:
        data = pandas_data_frame

    # data is in the form from te_clounts/scTE
    # rows = barcode
    # cols = genes

    print('Loaded Data Frame')
    if ensg_to_symbol: # Fix the table so it has gene names in the prime slot
        gene_ensg = data.columns
        gene_names = []

        for ensg in gene_ensg:
            if ensg not in ensg_to_symbol:
                gene_names.append(ensg)
                if '?' not in ensg: # Only print the warning if not a ? marked TE
                    print(f'Warning: {ensg} not found')
            else:
                gene_names.append(ensg_to_symbol[ensg])

    else: # Just use the inbuilt labels;
        gene_names = data.columns
        gene_ensg = gene_names

    if drop_fusions or drop_mir or drop_tes:
        todrop, gene_names, gene_ensg = _drop_fusions_mir(data,
            gene_names, gene_ensg, ensg_to_symbol,
            drop_fusions, drop_mir, drop_tes, drop_ribosomes)
        print('Dropped {} fusions/mirs/tes'.format(len(todrop)))

    layers = None
    if velocyte_data:
        _number_of_matched_genes = 0 # Number of genes matched in the velocyto and scTE/te_counts
        print('Velocyte data found, loading')
        spliced, unspliced, ambiguous, vel_genes, vel_barcodes = _load_velocyte_mtx(velocyte_data, load_ambiguous)

        print('Velocyte data fixing and matching barcodes genes')

        # First I need to get the rows (barcodes) in data
        vel_barcodes = {b: i for i, b in enumerate(list(vel_barcodes.index))}
        vel_genes = {g: i for i, g in enumerate(list(vel_genes.index))}

        barcode_indeces_to_keep = []
        for bc in list(data.index):
            # Possible to have a missing barcode?
            barcode_indeces_to_keep.append(vel_barcodes[bc])

        spliced = spliced[barcode_indeces_to_keep, :]
        unspliced = unspliced[barcode_indeces_to_keep, :]
        if load_ambiguous: ambiguous = ambiguous[barcode_indeces_to_keep, :]

        # stick a dummy 'empty' gene on the end I can use as a TE or missing gene
        np.vstack((spliced, np.zeros([len(barcode_indeces_to_keep), 1])))
        np.vstack((unspliced, np.zeros([len(barcode_indeces_to_keep), 1])))
        if load_ambiguous: np.vstack((ambiguous, np.zeros([len(barcode_indeces_to_keep), 1])))

        index_of_dummy_TE = spliced.shape[1]-1

        # Now the same for the columns/genes
        gene_indeces_to_keep = []
        need_to_fill_in = []
        for gene in gene_names:
            # Possible to have a missing barcode?
            if gene not in vel_genes:
                gene_indeces_to_keep.append(index_of_dummy_TE)
                continue

            _number_of_matched_genes += 1
            gene_indeces_to_keep.append(vel_genes[gene])

        print(f'Number of matched genes from scTE/te_counts and Velocyto: {_number_of_matched_genes}')

        spliced = spliced[:,gene_indeces_to_keep]
        unspliced = unspliced[:,gene_indeces_to_keep]
        if load_ambiguous: ambiguous = ambiguous[:,gene_indeces_to_keep]

        print('Done Velocyte')
        if load_ambiguous:
            layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous}
        else:
            layers = {'spliced': spliced, 'unspliced': unspliced}

    cells = data.index
    print('Sparsifying {}'.format(data.shape))
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')

    #print(len(cells), len(gene_names))
    #print(spliced.shape)
    #print(unspliced.shape)
    #print(ambiguous.shape)
    #print(data.shape)

    print('Building AnnData object')
    ad = AnnData(
        data,
        obs={'obs_names': cells},
        var={'var_names': gene_names, 'names': gene_ensg},
        layers=layers,
        dtype='float32')

    del data

    ad.var_names_make_unique()

    for k in obs_add:
        ad.obs[k] = obs_add[k]

    e = time.time()
    print('Done {:.1f}'.format(e-s))
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

def merge_barcode_umi_fastqs(barcode_fastq, umi_fastq, output_fastq, gzip=True):
    """
    **Purpose**
        Some FASTQ data (especially early version1 10x data) has the barcode and UMIs in a separate file.
        But STARsolo and newer cellranger (and other tools) require the CB and UMI to be merged into a
        single FASTQ.

    **Arguments**
        barcode_fastq (Required)
        umi_fastq (Required)
        output_fastq (Required)

        gzip (Optional, default=True)
            THe FASTQ files (and the output) are gziped.

    **Returns**
        THe lengths of the barcode and UMIs in a dictionary
    """
    assert barcode_fastq, 'You muse specify a barcode_fastq filename'
    assert umi_fastq, 'You muse specify a umi_fastq filename'
    assert output_fastq, 'You muse specify a output_fastq filename'
    assert os.path.exists(barcode_fastq), '{0} barcode_fastq not found'.format(barcode_fastq)
    assert os.path.exists(umi_fastq), '{0} barcode_fastq not found'.format(umi_fastq)

    print('Started {}'.format(output_fastq))

    barcode_len = 0
    umi_len = 0

    if gzip:
        output = gzipfile.open(output_fastq, 'wt')
    else:
        output = open(output_fastq, 'wt')

    for idx, read in enumerate(fastqPE(barcode_fastq, umi_fastq)):
        if (idx+1) % 100000 == 0:
            print('Done {:,}'.format(idx+1))

        output.write('{0}\n'.format(read[0]['name']))
        output.write('{0}{1}\n'.format(read[0]['seq'], read[1]['seq']))
        output.write('{0}\n'.format(read[0]['strand']))
        output.write('{0}{1}\n'.format(read[0]['qual'], read[1]['qual']))

    output.close()

    return {'barcode_len': barcode_len, 'umi_len': umi_len}

# From:
# https://github.com/alexdobin/STAR/issues/774#issuecomment-850477636
def buildAnnDataFromStarForscVelo(path):
    """
    Generate an anndata object from the STAR aligner output folder
    """
    path=path
    print('Load Read Counts')
    X = sc.read_mtx(os.path.join(path, 'Gene/raw/matrix.mtx'))

    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    X = X.X.transpose()

    print('Loading Spliced, Unspliced and Ambiguous matrices')
    mtxU = np.loadtxt(os.path.join(path, 'Velocyto/raw/unspliced.mtx'), skiprows=3, delimiter=' ')
    mtxS = np.loadtxt(os.path.join(path, 'Velocyto/raw/spliced.mtx'), skiprows=3, delimiter=' ')
    mtxA = np.loadtxt(os.path.join(path, 'Velocyto/raw/ambiguous.mtx'), skiprows=3, delimiter=' ')

    # Extract sparse matrix shape informations from the third row
    shapeU = np.loadtxt(os.path.join(path, 'Velocyto/raw/unspliced.mtx'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)
    shapeS = np.loadtxt(os.path.join(path, 'Velocyto/raw/spliced.mtx'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)
    shapeA = np.loadtxt(os.path.join(path, 'Velocyto/raw/ambiguous.mtx'), skiprows=2, max_rows=1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects

    spliced = sp.sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
    unspliced = sp.sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()
    ambiguous = sp.sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()

    print('Loading Genes and Identifiers')
    obs = pd.read_csv(os.path.join(path, 'Velocyto/raw/barcodes.tsv'), header=None, index_col=0)

    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None

    var = pd.read_csv(os.path.join(path, 'Velocyto/raw/features.tsv'), sep='\t',
        names=('gene_ids', 'feature_types'), index_col=1)

    print('Build AnnData object to be used with ScanPy and ScVelo')
    adata = AnnData(X=X, obs=obs, var=var,
        layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous})
    adata.var_names_make_unique()

    # Subset Cells based on STAR filtering
    selected_barcodes = pd.read_csv(os.path.join(path, 'Gene/filtered/barcodes.tsv'), header = None)
    adata = adata[selected_barcodes[0]]

    return adata

def smartseq_to_sparse(
    count_files:list,
    velocyto_looms:list = None,
    barcode_names:list = None,
    extra_obs:dict = None,
    extra_vars:dict = None,
    gzip:bool = True,
    test:int = None,
    ):
    '''

    **Purpose**
        Load count files into a seurat matrix

    **Arguments**
        count_files (Required)
            list of count matrices lo load. It assumes the filenames are the barcode names,
            unless 'barcode_names' is provided which will set the barcode names.

        velocyto_looms (Optional)

        gzip (Optional, default=True)
            count_files are gzipped

    **Returns**
        An AnnData object
    '''
    print(f'Starting smartseq_to_sparse() with {len(count_files)} to parse')

    if barcode_names: assert len(barcode_names) == len(count_files), 'count_files and barcode_names must be the same length'

    ###### init features and barcodes

    features_seen = []
    if barcode_names:
        barcodes_seen = barcode_names
    else:
        barcode_names = [os.path.split(f)[1].replace('.gz', '').replace('.tsv', '').replace('.csv', '') for f in count_files]

    data_dict = {}

    p = progressbar(len(count_files))
    for bcidx, f in enumerate(count_files):
        if test and bcidx > test:
            break

        # Load a table;
        features = []
        counts = []

        oh = open(f, 'rt') if not gzip else gzipfile.open(f, 'rt')
        for line in oh:
            line = line.strip().split('\t')

            # scTE/te_counts format;
            features.append(line[0])
            counts.append(int(line[1]))

        if not features_seen:
            features_seen = features # Assume first seen is the definitive version
            features_seen_set = set(features_seen)
        else:
            # Make sure features match
            assert len(features_seen_set) == len(set(features)), f'Mismatch in the length of the features_seen and the features in this sample {f}'
            assert features_seen_set == set(features), f'The features_seen and the features in this sample {f} do not match'

        bc = barcode_names[bcidx]
        data_dict[bc] = counts

        oh.close()
        p.update(bcidx)

    print('Loaded matrices')

    # put data_dict into a matrix

    if test:
        m = np.array(data_dict[barcode_names[0]])
        barcodes_used = [barcode_names[0],]
        for bc in barcode_names[1:]: # This could be done all in one go;
            if bc in data_dict: # in case test is being used.
                m = np.vstack((m, data_dict[bc]))
                barcodes_used.append(bc)
    else: # Much faster
        m = np.vstack(tuple(data_dict.values()))
        barcodes_used = barcode_names

    print('Build numpy matrix')

    # columns are barcodes
    # rows are features

    print()
    print(m)

    m = sp.sparse.csr_matrix(m)
    print('Sparsified')

    obs = {'obs_names': barcodes_used}
    if extra_obs:
        obs.update(extra_obs)

    var = {'var_names': features_seen}
    if extra_vars:
        if 'var_names' in extra_vars: # override if present;
            var = {}
            var['features_seen'] = features_seen
        var.update(extra_vars)

    ad = AnnData(X=m,
        obs=obs,
        var=var, #, 'names': gene_ensg},
        dtype='float')

    ad.var_names_make_unique()

    ###### Load  velocyto

    return ad

