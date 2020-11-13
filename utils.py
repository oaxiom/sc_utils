
import sys, os
import gzip as gzipfile
import scipy as sp
import pandas as pd
from anndata import AnnData

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


def sparsify(filename, obs_add, csv=True, drop_fusions=False, drop_mir=False):
    '''
    **Purpose**
        Convert a dense array in filename into a sparse array and return a
        AnnData object

    '''
    print('Started {0}'.format(filename))
    if csv:
        data = pd.read_csv(filename, index_col=0, header=0)
    else:
        data = pd.read_csv(filename, index_col=0, header=0, sep='\t')

    genes = data.columns

    todrop = []
    if drop_fusions:
        # drop genes with - in the form, but not -AS
        for i in genes:
            if '-' in i:
                if '-AS' in i:
                    continue
                if '-int' in i:
                    continue
                todrop.append(i)

        data.drop(todrop, axis=1)
        print('Dropped {} fusions'.format(len(todrop)))

    if drop_mir:
        todrop = [i for i in genes if i[0:3] == 'MIR']
        data.drop(todrop, axis=1)
        print('Dropped {} MIR'.format(len(todrop)))

    cells = data.index
    print('Sparsifying')
    data = sp.sparse.csr_matrix(data.to_numpy())
    data.astype('float32')

    '''
    oh = open('gene_names.{0}.tsv'.format(os.path.split(filename)[1]), 'w')
    for g in genes:
        oh.write('%s\n' % g)
    oh.close()
    '''

    print('Loaded')
    ad = AnnData(data, obs={'obs_names': cells}, var={'var_names': genes})
    del data

    for k in obs_add:
        ad.obs[k] = obs_add[k]

    print('Done')
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
        But STARsolo and newer cellreanger (and other tools) require the CB and UMI to be merged into a
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

    print('Started {0}'.format(output_fastq))

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
