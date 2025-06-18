import networkx as nx
import numpy as np
import pandas as pd

from anndata import AnnData

from scanpy import logging as logg
from scipy.sparse import csr_matrix

from dandelion.utilities._core import Dandelion
from dandelion.utilities._utilities import MuData, Tree, type_check


def transfer(
    adata: AnnData,
    dandelion: Dandelion,
    expanded_only: bool = False,
    neighbors_key: str | None = None,
    rna_key: str | None = None,
    vdj_key: str | None = None,
    clone_key: str | None = None,
    collapse_nodes: bool = False,
    overwrite: bool | list[str] | str | None = None,
) -> None:
    """
    Transfer data in `Dandelion` slots to `AnnData` object, updating the `.obs`, `.uns`, `.obsm` and `.obsp`slots.

    Parameters
    ----------
    adata : AnnData
        `AnnData` object.
    dandelion : Dandelion
        `Dandelion` object.
    expanded_only : bool, optional
        Whether or not to transfer the embedding with all cells with BCR (False) or only for expanded clones (True).
    neighbors_key : str | None, optional
        key for 'neighbors' slot in `.uns`.
    rna_key : str | None, optional
        prefix for stashed RNA connectivities and distances.
    vdj_key : str | None, optional
        prefix for stashed VDJ connectivities and distances.
    clone_key : str | None, optional
        column name of clone/clonotype ids. Only used for integration with scirpy.
    collapse_nodes : bool, optional
        Whether or not to transfer a cell x cell or clone x clone connectivity matrix into `.uns`. Only used for
        integration with scirpy.
    overwrite : bool | list[str] | str | None, optional
        Whether or not to overwrite existing anndata columns. Specifying a string indicating column name or
        list of column names will overwrite that specific column(s).
    """
    start = logg.info("Transferring network")
    # always overwrite with whatever columns are in dandelion's metadata:
    for x in dandelion.metadata.columns:
        if x not in adata.obs.columns:
            adata.obs[x] = pd.Series(dandelion.metadata[x])
        elif overwrite is True:
            adata.obs[x] = pd.Series(dandelion.metadata[x])
        if type_check(dandelion.metadata, x):
            adata.obs[x] = adata.obs[x].replace(np.nan, "No_contig")
        if adata.obs[x].dtype == "bool":
            adata.obs[x] = [str(x) for x in adata.obs[x]]

    if (overwrite is not None) and (overwrite is not True):
        if not type(overwrite) is list:
            overwrite = [overwrite]
        for ow in overwrite:
            adata.obs[ow] = pd.Series(dandelion.metadata[ow])
            if type_check(dandelion.metadata, ow):
                adata.obs[ow] = adata.obs[ow].replace(np.nan, "No_contig")

    if dandelion.graph is not None:
        if expanded_only:
            G = dandelion.graph[1]
        else:
            G = dandelion.graph[0]
        logg.info("converting matrices")
        distances = nx.to_pandas_adjacency(
            G, dtype=np.float32, weight="weight", nonedge=np.nan
        )
        df_distances = distances.reindex(
            index=adata.obs_names, columns=adata.obs_names
        )
        # convert to connectivity
        distances = distances.apply(lambda x: np.maximum(1e-45, 1 / np.exp(x)))
        df_connectivities = distances.reindex(
            index=adata.obs_names, columns=adata.obs_names
        )

        df_connectivities = df_connectivities.values
        df_connectivities[np.isnan(df_connectivities)] = 0
        df_distances = df_distances.values
        df_distances[np.isnan(df_distances)] = 0

        df_connectivities_ = csr_matrix(df_connectivities, dtype=np.float32)
        df_distances = csr_matrix(df_distances, dtype=np.float32)

        logg.info("Updating anndata slots")
        if neighbors_key is None:
            neighbors_key = "neighbors"
        if neighbors_key not in adata.uns:
            skip_stash = True
            adata.uns[neighbors_key] = {}

        rna_neighbors_key = "rna_" + neighbors_key
        vdj_neighbors_key = "vdj_" + neighbors_key
        if rna_neighbors_key not in adata.uns:
            adata.uns[rna_neighbors_key] = adata.uns[neighbors_key].copy()

        if rna_key is None:
            r_connectivities_key = "rna_connectivities"
            r_distances_key = "rna_distances"
        else:
            r_connectivities_key = rna_key + "_connectivitites"
            r_distances_key = rna_key + "_distances"

        if vdj_key is None:
            b_connectivities_key = "vdj_connectivities"
            b_distances_key = "vdj_distances"
        else:
            b_connectivities_key = vdj_key + "_connectivitites"
            b_distances_key = vdj_key + "_distances"

        # stash_rna_connectivities:
        if r_connectivities_key not in adata.obsp:
            if "skip_stash" not in locals():
                try:
                    adata.obsp[r_connectivities_key] = adata.obsp[
                        "connectivities"
                    ].copy()
                    adata.obsp[r_distances_key] = adata.obsp["distances"].copy()
                except:
                    adata.obsp[r_connectivities_key] = adata.uns[neighbors_key][
                        "connectivities"
                    ]
                    adata.obsp[r_distances_key] = adata.uns[neighbors_key][
                        "distances"
                    ]

        # always overwrite the bcr slots
        adata.obsp["connectivities"] = df_connectivities_.copy()
        adata.obsp["distances"] = df_distances.copy()
        adata.obsp[b_connectivities_key] = adata.obsp["connectivities"].copy()
        adata.obsp[b_distances_key] = adata.obsp["distances"].copy()

        # create the dictionary that will enable the use of scirpy's plotting.
        clonekey = clone_key if clone_key is not None else "clone_id"

        if not collapse_nodes:
            df_connectivities_[df_connectivities_.nonzero()] = 1
            cell_indices = {
                str(i): np.array([k])
                for i, k in zip(range(0, len(adata.obs_names)), adata.obs_names)
            }
        else:
            cell_indices = Tree()
            for x, y in adata.obs[clonekey].items():
                if y not in [
                    "",
                    "unassigned",
                    np.nan,
                    "NaN",
                    "NA",
                    "nan",
                    "None",
                    None,
                    "none",
                ]:
                    cell_indices[y][x].value = 1
            cell_indices = {
                str(x): np.array(list(r))
                for x, r in zip(
                    range(0, len(cell_indices)), cell_indices.values()
                )
            }
            df_connectivities_ = np.zeros(
                [len(cell_indices), len(cell_indices)]
            )
            np.fill_diagonal(df_connectivities_, 1)
            df_connectivities_ = csr_matrix(df_connectivities_)

        adata.uns[clonekey] = {
            # this is a symmetrical, pairwise, sparse distance matrix of clonotypes
            # the matrix is offset by 1, i.e. 0 = no connection, 1 = distance 0
            "distances": df_connectivities_,
            # '0' refers to the row/col index in the `distances` matrix
            # (numeric index, but needs to be strbecause of h5py)
            # np.array(["cell1", "cell2"]) points to the rows in `adata.obs`
            "cell_indices": cell_indices,
        }

    tmp = adata.obs.copy()
    if dandelion.graph is not None:
        if dandelion.layout is not None:
            if expanded_only:
                coord = pd.DataFrame.from_dict(
                    dandelion.layout[1], orient="index"
                )
            else:
                coord = pd.DataFrame.from_dict(
                    dandelion.layout[0], orient="index"
                )
            for x in coord.columns:
                tmp[x] = coord[x]

            X_vdj = np.array(tmp[[0, 1]], dtype=np.float32)
            adata.obsm["X_vdj"] = X_vdj

        logg.info(
            " finished",
            time=start,
            deep=(
                "updated `.obs` with `.metadata`\n"
                "added to `.uns['"
                + neighbors_key
                + "']` and `.uns['"
                + clonekey
                + "']`\n"
                "and `.obsp`\n"
                "   'distances', clonotype-weighted adjacency matrix\n"
                "   'connectivities', clonotype-weighted adjacency matrix"
            ),
        )
    else:
        logg.info(
            " finished", time=start, deep=("updated `.obs` with `.metadata`\n")
        )


tf = transfer  # alias for transfer function


def from_ak(airr: "Array") -> pd.DataFrame:
    """
    Convert an AIRR-formatted array to a pandas DataFrame.

    Parameters
    ----------
    airr : Array
        The AIRR-formatted array to be converted.

    Returns
    -------
    pd.DataFrame
        The converted pandas DataFrame.

    Raises
    ------
    KeyError
        If `sequence_id` not found in the data.
    """
    import awkward as ak

    df = ak.to_dataframe(airr)
    # check if 'sequence_id' column does not exist or if any value in 'sequence_id' is NaN
    if "sequence_id" not in df.columns or df["sequence_id"].isnull().any():
        df_reset = df.reset_index()

        # create a new 'sequence_id' column
        df_reset["sequence_id"] = df_reset.apply(
            lambda row: f"{row['cell_id']}_contig_{row['subentry'] + 1}", axis=1
        )

        # set 'entry' and 'subentry' back as the index
        df = df_reset.set_index(["entry", "subentry"])

    if "sequence_id" in df.columns:
        df.set_index("sequence_id", drop=False, inplace=True)
    if "cell_id" not in df.columns:
        df["cell_id"] = [c.split("_contig")[0] for c in df["sequence_id"]]

    return df


def to_ak(
    data: pd.DataFrame,
    **kwargs,
) -> tuple["Array", pd.DataFrame]:
    """
    Convert data from a DataFrame to an AnnData object with AIRR format.

    Parameters
    ----------
    data : pd.DataFrame
        The input DataFrame containing the data.
    **kwargs
        Additional keyword arguments passed to `scirpy.io.read_airr`.

    Returns
    -------
    tuple[Array, pd.DataFrame]
        A tuple containing the AIRR-formatted data as an ak.Array and the cell-level attributes as a pd.DataFrame.
    """

    try:
        import scirpy as ir
    except:
        raise ImportError("Please install scirpy to use this function.")

    adata = ir.io.read_airr(data, **kwargs)

    return adata.obsm["airr"], adata.obs


def _create_anndata(
    airr: "Array",
    obs: pd.DataFrame,
    adata: AnnData | None = None,
) -> AnnData:
    """
    Create an AnnData object with the given AIRR array and observation data.

    Parameters
    ----------
    airr : Array
        The AIRR array.
    obs : pd.DataFrame
        The observation data.
    adata : AnnData | None, optional
        An existing AnnData object to update. If None, a new AnnData object will be created.

    Returns
    -------
    AnnData
        The AnnData object with the AIRR array and observation data.
    """
    obsm = {"airr": airr}
    temp = AnnData(X=None, obs=obs, obsm=obsm)

    if adata is None:
        adata = temp
    else:
        cell_names = adata.obs_names.intersection(temp.obs_names)
        adata = adata[adata.obs_names.isin(cell_names)].copy()
        temp = temp[temp.obs_names.isin(cell_names)].copy()
        adata.obsm = dict() if adata.obsm is None else adata.obsm
        adata.obsm.update(temp.obsm)

    return adata


def _create_mudata(
    gex: AnnData,
    adata: AnnData,
    key: tuple[str, str] = ("gex", "airr"),
) -> MuData:
    """
    Create a MuData object from the given AnnData objects.

    Parameters
    ----------
    gex : AnnData
        The AnnData object containing gene expression data.
    adata : AnnData
        The AnnData object containing additional data.
    key : tuple[str, str], optional
        The keys to use for the gene expression and additional data in the MuData object. Defaults to ("gex", "airr").

    Returns
    -------
    MuData
        The created MuData object.

    Raises
    ------
    ImportError
        If the mudata package is not installed.
    """

    try:
        import mudata
    except ImportError:
        raise ImportError("Please install mudata. pip install mudata")
    if gex is not None:
        return mudata.MuData({key[0]: gex, key[1]: adata})
    return mudata.MuData({key[1]: adata})


def to_scirpy(
    data: Dandelion,
    transfer: bool = False,
    to_mudata: bool = True,
    gex_adata: AnnData | None = None,
    key: tuple[str, str] = ("gex", "airr"),
    **kwargs,
) -> AnnData | MuData:
    """
    Convert Dandelion data to scirpy-compatible format.

    Parameters
    ----------
    data : Dandelion
        The Dandelion object containing the data to be converted.
    transfer : bool, optional
        Whether to transfer additional information from Dandelion to the converted data. Defaults to False.
    to_mudata : bool, optional
        Whether to convert the data to MuData format instead of AnnData. Defaults to True.
        If converting to AnnData, it will assert that the same cell_ids and .obs_names are present in the `gex_adata` provided.
    gex_adata : AnnData, optional
        An existing AnnData object to be used as the base for the converted data if provided.
    key : tuple[str, str], optional
        A tuple specifying the keys for the 'gex' and 'airr' fields in the converted data. Defaults to ("gex", "airr").
    **kwargs
        Additional keyword arguments passed to `scirpy.io.read_airr`.

    Returns
    -------
    AnnData | MuData
        The converted data in either AnnData or MuData format.
    """

    if "umi_count" not in data.data and "duplicate_count" in data.data:
        data.data["umi_count"] = data.data["duplicate_count"]
    for h in [
        "sequence",
        "rev_comp",
        "sequence_alignment",
        "germline_alignment",
        "v_cigar",
        "d_cigar",
        "j_cigar",
    ]:
        if h not in data.data:
            data.data[h] = None

    airr, obs = to_ak(data.data, **kwargs)
    if to_mudata:
        adata = _create_anndata(
            airr,
            obs,
        )
        if transfer:
            tf(adata, data)  # need to make a version that is not so verbose?

        return _create_mudata(gex_adata, adata, key)
    else:
        adata = _create_anndata(airr, obs, gex_adata)

        if transfer:
            tf(adata, data)
        return adata


def from_scirpy(data: AnnData | MuData) -> Dandelion:
    """
    Convert data from scirpy format to Dandelion format.

    Parameters
    ----------
    data : AnnData | MuData
        The input data in scirpy format.

    Returns
    -------
    Dandelion
        The converted data in Dandelion format.
    """

    if not isinstance(data, AnnData):
        data = data.mod["airr"]
    data = data.copy()
    data.obsm["airr"]["cell_id"] = data.obs.index
    data = from_ak(data.obsm["airr"])

    return Dandelion(data)
