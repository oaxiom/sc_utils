
from .utils import (sparsify,
                    export_dense,
                    merge_barcode_umi_fastqs,
                    cmap_grey_red,
                    buildAnnDataFromStarForscVelo,
                    smartseq_to_sparse,
                    cell_type_prop_bar)

from pkgutil import iter_modules
available_modules = list((name for loader, name, ispkg in iter_modules()))

if 'cvxpy' in available_modules:
    from .scran import (compute_sum_factors)
