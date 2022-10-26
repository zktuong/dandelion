#!/usr/bin/env python
"""test trajectory"""
import dandelion as ddl
import scanpy as sc
import numpy as np
import urllib.request
import pandas as pd
import sys
import pytest
from unittest.mock import patch


@pytest.mark.skipif(
    (sys.platform == "darwin") & (sys.version_info.minor < 8),
    reason="macos CI stalls.",
)
@patch("matplotlib.pyplot.show")
def test_trajectory(mock_show):
    """test_workflow"""
    import milopy.core as milo
    import palantir

    file = "demo-pseudobulk.h5ad"
    fname = "ftp://ftp.sanger.ac.uk/pub/users/kp9/" + file
    urllib.request.urlretrieve(fname, file)
    adata = sc.read(file)
    adata = ddl.tl.setup_vdj_pseudobulk(adata)
    sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=50)
    milo.make_nhoods(adata)
    sc.tl.umap(adata)
    pb_adata = ddl.tl.vdj_pseudobulk(
        adata, pbs=adata.obsm["nhoods"], obs_to_take="anno_lvl_2_final_clean"
    )
    sc.tl.pca(pb_adata)
    sc.pl.pca(pb_adata, color="anno_lvl_2_final_clean")
    # palantir business
    rootcell = np.argmax(pb_adata.obsm["X_pca"][:, 0])
    terminal_states = pd.Series(
        ["CD8+T", "CD4+T"],
        index=pb_adata.obs_names[
            [
                np.argmax(pb_adata.obsm["X_pca"][:, 1]),
                np.argmin(pb_adata.obsm["X_pca"][:, 1]),
            ]
        ],
    )
    # Run diffusion maps
    pca_projections = pd.DataFrame(
        pb_adata.obsm["X_pca"], index=pb_adata.obs_names
    )
    dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=5)
    ms_data = palantir.utils.determine_multiscale_space(dm_res)
    pr_res = palantir.core.run_palantir(
        ms_data,
        pb_adata.obs_names[rootcell],
        num_waypoints=500,
        terminal_states=terminal_states.index,
    )
    pr_res.branch_probs.columns = terminal_states[pr_res.branch_probs.columns]
    pb_adata = ddl.tl.pseudotime_transfer(pb_adata, pr_res)
    bdata = ddl.tl.project_pseudotime_to_cell(
        adata, pb_adata, terminal_states.values
    )
