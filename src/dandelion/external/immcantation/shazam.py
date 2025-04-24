import os
import re
import functools
import warnings

import pandas as pd
import numpy as np

from scanpy import logging as logg
from typing import Literal

from plotnine import (
    ggplot,
    options,
    aes,
    xlab,
    ylab,
    facet_wrap,
    theme,
    annotate,
    theme_bw,
    geom_histogram,
    geom_vline,
    save_as_pdf_pages,
)

from dandelion.utilities._core import Dandelion
from dandelion.utilities._utilities import load_data, write_airr, sanitize_data


def quantify_mutations(
    data: Dandelion | str,
    split_locus: bool = False,
    sequence_column: str | None = None,
    germline_column: str | None = None,
    region_definition: str | None = None,
    mutation_definition: str | None = None,
    frequency: bool = False,
    combine: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    Run basic mutation load analysis.

    Implemented in `shazam` https://shazam.readthedocs.io/en/stable/vignettes/Mutation-Vignette.

    Parameters
    ----------
    data : Dandelion | str
        `Dandelion` object, file path to AIRR file.
    split_locus : bool, optional
        whether to return the results for heavy chain and light chain separately.
    sequence_column : str | None, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/observedMutations
    germline_column : str | None, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/observedMutations
    region_definition : str | None, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/IMGT_SCHEMES/
    mutation_definition : str | None, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/MUTATION_SCHEMES/
    frequency : bool, optional
        whether to return the results a frequency or counts.
    combine : bool, optional
        whether to return the results for replacement and silent mutations separately.
    **kwargs
        passed to shazam::observedMutations.

    Returns
    -------
    pd.DataFrame
        pandas DataFrame holding mutation information.
    """
    start = logg.info("Quantifying mutations")
    try:
        import rpy2
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import NULL
        from rpy2.robjects import pandas2ri
    except:
        raise (
            ImportError(
                "Unable to initialise R instance. Please run this separately through R with shazam's tutorials."
            )
        )

    sh = importr("shazam")
    base = importr("base")
    if isinstance(data, Dandelion):
        dat = load_data(data.data)
    else:
        dat = load_data(data)

    pandas2ri.activate()
    warnings.filterwarnings("ignore")

    dat = sanitize_data(dat)

    if "ambiguous" in dat:
        dat_ = dat[dat["ambiguous"] == "F"].copy()
    else:
        dat_ = dat.copy()

    if sequence_column is None:
        seq_ = "sequence_alignment"
    else:
        seq_ = sequence_column

    if germline_column is None:
        germline_ = "germline_alignment_d_mask"
    else:
        germline_ = germline_column

    if region_definition is None:
        reg_d = NULL
    else:
        reg_d = base.get(region_definition)

    if mutation_definition is None:
        mut_d = NULL
    else:
        mut_d = base.get(mutation_definition)

    if split_locus is False:
        dat_ = dat_.where(dat_.isna(), dat_.astype(str))
        try:
            dat_r = pandas2ri.py2rpy(dat_)
        except:
            dat_ = dat_.astype(str)
            dat_r = pandas2ri.py2rpy(dat_)

        results = sh.observedMutations(
            dat_r,
            sequenceColumn=seq_,
            germlineColumn=germline_,
            regionDefinition=reg_d,
            mutationDefinition=mut_d,
            frequency=frequency,
            combine=combine,
            **kwargs,
        )
        if rpy2.__version__ >= "3.4.5":
            from rpy2.robjects.conversion import localconverter

            with localconverter(
                rpy2.robjects.default_converter + pandas2ri.converter
            ):
                pd_df = rpy2.robjects.conversion.rpy2py(results)
        else:
            # pd_df = pandas2ri.rpy2py_data frame(results)
            pd_df = results.copy()
    else:
        dat_h = dat_[dat_["locus"] == "IGH"]
        dat_l = dat_[dat_["locus"].isin(["IGK", "IGL"])]

        dat_h = dat_h.where(dat_h.isna(), dat_h.astype(str))
        try:
            dat_h_r = pandas2ri.py2rpy(dat_h)
        except:
            dat_h = dat_h.astype(str)
            dat_h_r = pandas2ri.py2rpy(dat_h)

        dat_l = dat_l.where(dat_l.isna(), dat_l.astype(str))
        try:
            dat_l_r = pandas2ri.py2rpy(dat_l)
        except:
            dat_l = dat_l.astype(str)
            dat_l_r = pandas2ri.py2rpy(dat_l)

        results_h = sh.observedMutations(
            dat_h_r,
            sequenceColumn=seq_,
            germlineColumn=germline_,
            regionDefinition=reg_d,
            mutationDefinition=mut_d,
            frequency=frequency,
            combine=combine,
            **kwargs,
        )
        results_l = sh.observedMutations(
            dat_l_r,
            sequenceColumn=seq_,
            germlineColumn=germline_,
            regionDefinition=reg_d,
            mutationDefinition=mut_d,
            frequency=frequency,
            combine=combine,
            **kwargs,
        )
        if rpy2.__version__ >= "3.4.5":
            from rpy2.robjects.conversion import localconverter

            with localconverter(
                rpy2.robjects.default_converter + pandas2ri.converter
            ):
                results_h = rpy2.robjects.conversion.rpy2py(results_h)
                results_l = rpy2.robjects.conversion.rpy2py(results_l)
        pd_df = pd.concat([results_h, results_l])

    pd_df.set_index("sequence_id", inplace=True, drop=False)
    # this doesn't actually catch overwritten columns
    cols_to_return = pd_df.columns.difference(dat.columns)
    if len(cols_to_return) < 1:
        cols_to_return = list(
            filter(re.compile("mu_.*").match, [c for c in pd_df.columns])
        )
    else:
        cols_to_return = cols_to_return

    res = {}
    if isinstance(data, Dandelion):
        for x in cols_to_return:
            res[x] = list(pd_df[x])
            # TODO: str will make it work for the back and forth conversion with rpy2. but maybe can use a better option
            dat_[x] = [str(r) for r in res[x]]
            data.data[x] = pd.Series(dat_[x])
        if split_locus is False:
            metadata_ = data.data[["cell_id"] + list(cols_to_return)]
        else:
            metadata_ = data.data[["locus", "cell_id"] + list(cols_to_return)]

        for x in cols_to_return:
            metadata_[x] = metadata_[x].astype(float)

        if split_locus is False:
            metadata_ = metadata_.groupby("cell_id").sum()
        else:
            metadata_ = metadata_.groupby(["locus", "cell_id"]).sum()
            metadatas = []
            for x in list(set(data.data["locus"])):
                tmp = metadata_.iloc[
                    metadata_.index.isin([x], level="locus"), :
                ]
                tmp.index = tmp.index.droplevel()
                tmp.columns = [c + "_" + str(x) for c in tmp.columns]
                metadatas.append(tmp)
            metadata_ = functools.reduce(
                lambda x, y: pd.merge(
                    x, y, left_index=True, right_index=True, how="outer"
                ),
                metadatas,
            )

        metadata_.index.name = None
        data.data = sanitize_data(data.data)
        if data.metadata is None:
            data.metadata = metadata_
        else:
            for x in metadata_.columns:
                data.metadata[x] = pd.Series(metadata_[x])
        logg.info(
            " finished",
            time=start,
            deep=(
                "Updated Dandelion object: \n"
                "   'data', contig-indexed AIRR table\n"
                "   'metadata', cell-indexed observations table\n"
            ),
        )
    else:
        for x in cols_to_return:
            res[x] = list(pd_df[x])
            # TODO: str will make it work for the back and forth conversion with rpy2. but maybe can use a better option
            dat[x] = [str(r) for r in res[x]]
        # dat = sanitize_data(dat)
        if isinstance(data, pd.DataFrame):
            logg.info(" finished", time=start, deep=("Returning DataFrame\n"))
            return dat
        elif os.path.isfile(data):
            logg.info(
                " finished",
                time=start,
                deep=("saving DataFrame at {}\n".format(str(data))),
            )
            write_airr(dat, data)


def calculate_threshold(
    data: Dandelion | pd.DataFrame | str,
    mode: Literal["single-cell", "heavy"] = "single-cell",
    manual_threshold: float | None = None,
    VJthenLen: bool = False,
    onlyHeavy: bool = False,
    model: (
        Literal[
            "ham",
            "aa",
            "hh_s1f",
            "hh_s5f",
            "mk_rs1nf",
            "hs1f_compat",
            "m1n_compat",
        ]
        | None
    ) = None,
    normalize_method: Literal["len"] | None = None,
    threshold_method: Literal["gmm", "density"] | None = None,
    edge: float | None = None,
    cross: list[float] | None = None,
    subsample: int | None = None,
    threshold_model: (
        Literal["norm-norm", "norm-gamma", "gamma-norm", "gamma-gamma"] | None
    ) = None,
    cutoff: Literal["optimal", "intersect", "user"] | None = None,
    sensitivity: float | None = None,
    specificity: float | None = None,
    plot: bool = True,
    plot_group: str | None = None,
    figsize: tuple[float, float] = (4.5, 2.5),
    save_plot: str | None = None,
    ncpu: int = 1,
    **kwargs,
) -> Dandelion:
    """
    Calculating nearest neighbor distances for tuning clonal assignment with `shazam`.

    https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/

    Runs the following:

    distToNearest
        Get non-zero distance of every heavy chain (IGH) sequence (as defined by sequenceColumn) to its nearest sequence
        in a partition of heavy chains sharing the same V gene, J gene, and junction length (VJL), or in a partition of
        single cells with heavy chains sharing the same heavy chain VJL combination, or of single cells with heavy and
        light chains sharing the same heavy chain VJL and light chain VJL combinations.
    findThreshold
        automtically determines an optimal threshold for clonal assignment of Ig sequences using a vector of nearest
        neighbor distances. It provides two alternative methods using either a Gamma/Gaussian Mixture Model fit
        (threshold_method="gmm") or kernel density fit (threshold_method="density").

    Parameters
    ----------
    data : Dandelion | pd.DataFrame | str
        input `Danelion`, AIRR data as pandas DataFrame or path to file.
    mode : Literal["single-cell", "heavy"], optional
        accepts one of "heavy" or "single-cell".
        Refer to https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette.
    manual_threshold : float | None, optional
        value to manually plot in histogram.
    VJthenLen : bool, optional
        logical value specifying whether to perform partitioning as a 2-stage process.
        If True, partitions are made first based on V and J gene, and then further split
        based on junction lengths corresponding to sequenceColumn.
        If False, perform partition as a 1-stage process during which V gene, J gene, and junction length
        are used to create partitions simultaneously.
        Defaults to False.
    onlyHeavy : bool, optional
        use only the IGH (BCR) or TRB/TRD (TCR) sequences for grouping. Only applicable to single-cell mode.
        See groupGenes for further details.
    model : Literal["ham", "aa", "hh_s1f", "hh_s5f", "mk_rs1nf", "hs1f_compat", "m1n_compat", ] | None, optional
        underlying SHM model, which must be one of "ham","aa","hh_s1f","hh_s5f","mk_rs1nf","hs1f_compat","m1n_compat".
    normalize_method : Literal["len"] | None, optional
        method of normalization. The default is "len", which divides the distance by the length of the sequence group.
        If "none" then no normalization if performed.
    threshold_method : Literal["gmm", "density"] | None, optional
        string defining the method to use for determining the optimal threshold. One of "gmm" or "density".
    edge : float | None, optional
        upper range as a fraction of the data density to rule initialization of Gaussian fit parameters.
        Default value is 0.9 (or 90). Applies only when threshold_method="density".
    cross : list[float] | None, optional
        supplementary nearest neighbor distance vector output from distToNearest for initialization of the Gaussian fit
        parameters. Applies only when method="gmm".
    subsample : int | None, optional
        maximum number of distances to subsample to before threshold detection.
    threshold_model : Literal["norm-norm", "norm-gamma", "gamma-norm", "gamma-gamma"] | None, optional
        allows the user to choose among four possible combinations of fitting curves: "norm-norm", "norm-gamma",
        "gamma-norm", and "gamma-gamma". Applies only when method="gmm".
    cutoff : Literal["optimal", "intersect", "user"] | None, optional
        method to use for threshold selection: the optimal threshold "optimal", the intersection point of the two fitted
        curves "intersect", or a value defined by user for one of the sensitivity or specificity "user". Applies only
        when method="gmm".
    sensitivity : float | None, optional
        sensitivity required. Applies only when method="gmm" and cutoff="user".
    specificity : float | None, optional
        specificity required. Applies only when method="gmm" and cutoff="user".
    plot : bool, optional
        whether or not to return plot.
    plot_group : str | None, optional
        determines the fill color and facets.
    figsize : tuple[float, float], optional
        size of plot.
    save_plot : str | None, optional
        if specified, plot will be save with this path.
    ncpu : int, optional
        number of cpus to run `distToNearest`. defaults to 1.
    **kwargs
        passed to shazam's `distToNearest <https://shazam.readthedocs.io/en/stable/topics/distToNearest/>`__.

    Returns
    -------
    Dandelion
        Dandelion object with `.threshold` slot filled.

    Raises
    ------
    ValueError
        if automatic thresholding failed.
    """
    start = logg.info("Calculating threshold")
    try:
        import rpy2
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import NULL
        from rpy2.robjects import pandas2ri, FloatVector
    except:
        raise (
            ImportError(
                "Unable to initialise R instance. Please run this separately through R with shazam's tutorials."
            )
        )
    if isinstance(data, Dandelion):
        dat = load_data(data.data)
    elif isinstance(data, pd.DataFrame) or os.path.isfile(str(data)):
        dat = load_data(data)
        warnings.filterwarnings("ignore")

    sh = importr("shazam")
    pandas2ri.activate()
    if "v_call_genotyped" in dat.columns:
        v_call = "v_call_genotyped"
    else:
        v_call = "v_call"
    if model is None:
        model_ = "ham"
    else:
        model_ = model
    if normalize_method is None:
        norm_ = "len"
    else:
        norm_ = normalize_method
    if threshold_method is None:
        threshold_method_ = "density"
    else:
        threshold_method_ = threshold_method
    if subsample is None:
        subsample_ = NULL
    else:
        subsample_ = subsample

    if mode == "heavy":
        dat_h = dat[dat["locus"].isin(["IGH", "TRB", "TRD"])].copy()
        try:
            dat_h_r = pandas2ri.py2rpy(dat_h)
        except:
            dat_h = dat_h.astype(str)
            dat_h_r = pandas2ri.py2rpy(dat_h)

        dist_ham = sh.distToNearest(
            dat_h_r, vCallColumn=v_call, model=model_, normalize=norm_, **kwargs
        )
    elif mode == "single-cell":
        try:
            dat_r = pandas2ri.py2rpy(dat)
        except:
            dat = dat.astype(str)
            dat_r = pandas2ri.py2rpy(dat)
        try:
            dist_ham = sh.distToNearest(
                dat_r,
                cellIdColumn="cell_id",
                locusColumn="locus",
                VJthenLen=VJthenLen,
                vCallColumn=v_call,
                onlyHeavy=onlyHeavy,
                normalize=norm_,
                model=model_,
                nproc=ncpu,
                **kwargs,
            )
        except:
            logg.info(
                "Rerun this after filtering. For now, switching to heavy mode."
            )
            dat_h = dat[dat["locus"].isin(["IGH", "TRB", "TRD"])].copy()
            try:
                dat_h_r = pandas2ri.py2rpy(dat_h)
            except:
                dat_h = dat_h.astype(str)
                dat_h_r = pandas2ri.py2rpy(dat_h)

            dist_ham = sh.distToNearest(
                dat_h_r,
                vCallColumn=v_call,
                model=model_,
                normalize=norm_,
                nproc=ncpu,
                **kwargs,
            )
    if rpy2.__version__ >= "3.4.5":
        from rpy2.robjects.conversion import localconverter

        with localconverter(
            rpy2.robjects.default_converter + pandas2ri.converter
        ):
            dist_ham = rpy2.robjects.conversion.rpy2py(dist_ham)
    # Find threshold using density method
    dist = np.array(dist_ham["dist_nearest"])
    if manual_threshold is None:
        if threshold_method_ == "density":
            if edge is None:
                edge_ = 0.9
            else:
                edge_ = edge
            dist_threshold = sh.findThreshold(
                FloatVector(dist[~np.isnan(dist)]),
                method=threshold_method_,
                subsample=subsample_,
                edge=edge_,
            )
            threshold = np.array(dist_threshold.slots["threshold"])[0]
            if np.isnan(threshold):
                logg.info(
                    "      Threshold method 'density' did not return with any values. Switching to method = 'gmm'."
                )
                threshold_method_ = "gmm"
                if threshold_model is None:
                    threshold_model_ = "gamma-gamma"
                else:
                    threshold_model_ = threshold_model
                if cross is None:
                    cross_ = NULL
                else:
                    cross_ = cross
                if cutoff is None:
                    cutoff_ = "optimal"
                else:
                    cutoff_ = cutoff
                if sensitivity is None:
                    sen_ = NULL
                else:
                    sen_ = sensitivity
                if specificity is None:
                    spc_ = NULL
                else:
                    spc_ = specificity
                dist_threshold = sh.findThreshold(
                    FloatVector(dist[~np.isnan(dist)]),
                    method=threshold_method_,
                    model=threshold_model_,
                    cross=cross_,
                    subsample=subsample_,
                    cutoff=cutoff_,
                    sen=sen_,
                    spc=spc_,
                )
                if rpy2.__version__ >= "3.4.5":
                    from rpy2.robjects.conversion import localconverter

                    with localconverter(
                        rpy2.robjects.default_converter + pandas2ri.converter
                    ):
                        dist_threshold = rpy2.robjects.conversion.rpy2py(
                            dist_threshold
                        )

                threshold = np.array(dist_threshold.slots["threshold"])[0]
        else:
            if threshold_model is None:
                threshold_model_ = "gamma-gamma"
            else:
                threshold_model_ = threshold_model
            if cross is None:
                cross_ = NULL
            else:
                cross_ = cross
            if cutoff is None:
                cutoff_ = "optimal"
            else:
                cutoff_ = cutoff
            if sensitivity is None:
                sen_ = NULL
            else:
                sen_ = sensitivity
            if specificity is None:
                spc_ = NULL
            else:
                spc_ = specificity
            dist_threshold = sh.findThreshold(
                FloatVector(dist[~np.isnan(dist)]),
                method=threshold_method_,
                model=threshold_model_,
                cross=cross_,
                subsample=subsample_,
                cutoff=cutoff_,
                sen=sen_,
                spc=spc_,
            )
            if rpy2.__version__ >= "3.4.5":
                from rpy2.robjects.conversion import localconverter

                with localconverter(
                    rpy2.robjects.default_converter + pandas2ri.converter
                ):
                    dist_threshold = rpy2.robjects.conversion.rpy2py(
                        dist_threshold
                    )
            threshold = np.array(dist_threshold.slots["threshold"])[0]
        if np.isnan(threshold):
            raise ValueError(
                "Automatic thresholding failed. Please visually inspect the resulting distribution fits"
                + " and choose a threshold value manually."
            )
        # dist_ham = pandas2ri.rpy2py_data frame(dist_ham)
        tr = threshold
    else:
        tr = manual_threshold

    if plot:
        options.figure_size = figsize
        if plot_group is None:
            plot_group = "sample_id"
        else:
            plot_group = plot_group

        p = (
            ggplot(dist_ham, aes("dist_nearest", fill=str(plot_group)))
            + theme_bw()
            + xlab("Grouped Hamming distance")
            + ylab("Count")
            + geom_histogram(binwidth=0.01)
            + geom_vline(
                xintercept=tr, linetype="dashed", color="blue", size=0.5
            )
            + annotate(
                "text",
                x=tr + 0.02,
                y=10,
                label="Threshold:\n" + str(np.around(tr, decimals=2)),
                size=8,
                color="Blue",
            )
            + facet_wrap("~" + str(plot_group), scales="free_y")
            + theme(legend_position="none")
        )
        if save_plot is not None:
            save_as_pdf_pages([p], filename=save_plot, verbose=False)
        p.show()
    else:
        logg.info(
            "Automatic Threshold : "
            + str(np.around(threshold, decimals=2))
            + "\n method = "
            + str(threshold_method_)
        )
    if isinstance(data, Dandelion):
        data.threshold = tr
        logg.info(
            " finished",
            time=start,
            deep=(
                "Updated Dandelion object: \n"
                "   'threshold', threshold value for tuning clonal assignment\n"
            ),
        )
    else:
        output = Dandelion(dat)
        output.threshold = tr
        return output
