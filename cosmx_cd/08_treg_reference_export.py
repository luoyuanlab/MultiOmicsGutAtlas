#!/usr/bin/env python
"""CD scRNA-seq Tregs: export of the label-transfer reference (raw counts and pseudotime bins).

The Treg pseudotime object is an AnnData file (cd_scrna/12_treg_pseudotime.py); Seurat label transfer runs in R
(09_treg_label_transfer.R). This step writes the reference as the two CSV files that 09 reads.

Panels made
  none (reference of the label transfer behind Fig 5, Supp 9 and the CITE-seq Treg labels)

Inputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/treg_pseudotime.h5ad (cd_scrna/12_treg_pseudotime.py)

Outputs
  /share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/objects/label_transfer/ (treg_raw_counts_rowfix.csv, treg_pt_bin_rowfix.csv)

Run order of sections
  treg_reference_csv

How to run
  python <this file>                          run every section in the order below
  RUN_SECTIONS=<id>[,<id>...] python <file>   run only the listed sections (still in file order)
  Each section runs one analysis step and reads its inputs from disk, so a single section can be run once its
  upstream outputs exist. Command-line arguments and environment variables used by a section are listed in its
  own docstring.
  Paths are the lab's absolute paths; edit them to match the local data location.
  Environment: see the README section "Environment" (OPENBLAS_CORETYPE=SkylakeX on Intel Sapphire Rapids nodes).
"""
import os as _os
_THIS_DIR = _os.path.dirname(_os.path.abspath(__file__))   # repository folder of this file (used for utils/ imports)
_RUN_SECTIONS = [_s.strip() for _s in _os.environ.get('RUN_SECTIONS', '').split(',') if _s.strip()]


def _run(section):
    """True if `section` should run: all sections run in file order unless RUN_SECTIONS=id1,id2,... is set."""
    if _RUN_SECTIONS and section not in _RUN_SECTIONS:
        return False
    print(f'===== section {section} =====', flush=True)
    return True


# ======================================================================================================================
# Section treg_reference_csv: label-transfer reference (scRNA-seq Treg raw counts and pseudotime bins) as CSV
# Panels: label-transfer reference (scRNA-seq Tregs with pseudotime bins)
# ======================================================================================================================
if _run('treg_reference_csv'):
    __file__ = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix/code/label_transfer_fig4/treg_reference_csv.py'   # working folder of this step (relative paths resolve from here)
    """Writes the label-transfer reference read by 09_treg_label_transfer.R (section treg_reference_build):
  - treg_raw_counts_rowfix.csv : cells x genes raw counts (layer 'counts'), first column = cell barcode
  - treg_pt_bin_rowfix.csv     : index = cell barcode, one column pt_bin (pseudotime bin 0-4)
Source: objects/treg_pseudotime.h5ad (2,941 Tregs). Prints the detection rate of FOXP3, IL2RA, CTLA4 and NKG7 in the reference.
Run with OPENBLAS_CORETYPE=SkylakeX /data/ydn4687/geneformer_env/bin/python."""
    import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
    D = '/share/fsmresfiles/UC/revision_jenny/cd_treg_rebuild_rowfix'
    OUT = f'{D}/objects/label_transfer'

    t = ad.read_h5ad(f'{D}/objects/treg_pseudotime.h5ad')
    X = t.layers['counts']
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    assert np.allclose(X, np.round(X)), 'counts layer is not integer'
    X = np.round(X).astype(np.int32)
    print('cells x genes:', X.shape)

    # marker detection in the reference Tregs
    for g in ['FOXP3', 'IL2RA', 'CTLA4', 'NKG7']:
        print(f'  detection {g}: {(X[:, t.var_names.get_loc(g)] > 0).mean():.3f}')

    pd.DataFrame(X, index=t.obs_names, columns=t.var_names).to_csv(f'{OUT}/treg_raw_counts_rowfix.csv')
    meta = pd.DataFrame({'pt_bin': t.obs['pt_bin'].astype(int).values}, index=t.obs_names)
    meta.to_csv(f'{OUT}/treg_pt_bin_rowfix.csv')
    print(meta['pt_bin'].value_counts().sort_index().to_string())
    print('wrote', f'{OUT}/treg_raw_counts_rowfix.csv', f'{OUT}/treg_pt_bin_rowfix.csv')
