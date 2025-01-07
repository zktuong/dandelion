#!/usr/bin/env python
import bz2
import gzip
import h5py
import json
import os
import pickle
import re
import shutil

import _pickle as cPickle
import networkx as nx
import numpy as np
import pandas as pd

from anndata import AnnData
from collections import defaultdict, OrderedDict
from pathlib import Path
from scanpy import logging as logg
from typing import Literal

from dandelion.tools._tools import transfer as tf
from dandelion.utilities._core import *
from dandelion.utilities._utilities import *

pickle.HIGHEST_PROTOCOL = 4


AIRR = [
    "cell_id",
    "sequence_id",
    "sequence",
    "sequence_aa",
    "productive",
    "complete_vdj",
    "vj_in_frame",
    "locus",
    "v_call",
    "d_call",
    "j_call",
    "c_call",
    "junction",
    "junction_aa",
    "consensus_count",
    "umi_count",
    "cdr3_start",
    "cdr3_end",
    "sequence_length_10x",
    "high_confidence_10x",
    "is_cell_10x",
    "fwr1_aa",
    "fwr1",
    "cdr1_aa",
    "cdr1",
    "fwr2_aa",
    "fwr2",
    "cdr2_aa",
    "cdr2",
    "fwr3_aa",
    "fwr3",
    "fwr4_aa",
    "fwr4",
    "clone_id",
    "raw_consensus_id_10x",
    "exact_subclonotype_id_10x",
]
CELLRANGER = [
    "barcode",
    "contig_id",
    "sequence",
    "aa_sequence",
    "productive",
    "full_length",
    "frame",
    "chain",
    "v_gene",
    "d_gene",
    "j_gene",
    "c_gene",
    "cdr3_nt",
    "cdr3",
    "reads",
    "umis",
    "cdr3_start",
    "cdr3_stop",
    "length",
    "high_confidence",
    "is_cell",
    "fwr1",
    "fwr1_nt",
    "cdr1",
    "cdr1_nt",
    "fwr2",
    "fwr2_nt",
    "cdr2",
    "cdr2_nt",
    "fwr3",
    "fwr3_nt",
    "fwr4",
    "fwr4_nt",
    "raw_clonotype_id",
    "raw_consensus_id",
    "exact_subclonotype_id",
]


def fasta_iterator(fh: str) -> tuple[str, str]:
    """Read in a fasta file as an iterator."""
    while True:
        line = fh.readline()
        if line.startswith(">"):
            break
    while True:
        header = line[1:-1].rstrip()
        sequence = fh.readline().rstrip()
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith(">"):
                break
            sequence += line.rstrip()
        yield (header, sequence)
        if not line:
            return


def read_pkl(filename: str = "dandelion_data.pkl.pbz2") -> Dandelion:
    """
    Read in and returns a `Dandelion` class saved using pickle format.

    Parameters
    ----------
    filename : str, optional
        path to `.pkl` file. Depending on the extension, it will try to unzip accordingly.

    Returns
    -------
    Dandelion
        saved `Dandelion` object in pickle format.
    """
    if isBZIP(str(filename)):
        data = bz2.BZ2File(filename, "rb")
        data = cPickle.load(data)
    elif isGZIP(str(filename)):
        data = gzip.open(filename, "rb")
        data = cPickle.load(data)
    else:
        with open(filename, "rb") as f:
            data = cPickle.load(f)
    return data


def read_h5ddl(filename: Path | str = "dandelion_data.h5ddl") -> Dandelion:
    """
    Read in and returns a `Dandelion` class from .h5ddl format.

    Parameters
    ----------
    filename : Path | str, optional
        path to `.h5ddl` file

    Returns
    -------
    Dandelion
        `Dandelion` object.

    Raises
    ------
    AttributeError
        if `data` not found in `.h5ddl` file.
    """
    data = load_data(_read_h5_group(filename, group="data"))
    metadata = _read_h5_group(filename, group="metadata")
    try:
        metadata_names = _read_h5_group(filename, group="metadata_names")
        metadata.index = metadata_names
    except KeyError:  # pragma: no cover
        pass

    try:
        g_0 = _read_h5_csr_matrix(filename, group="graph/graph_0")
        g_1 = _read_h5_csr_matrix(filename, group="graph/graph_1")
        graph0 = _create_graph(g_0, adjust_adjacency=1, fillna=0)
        graph1 = _create_graph(g_1, adjust_adjacency=1, fillna=0)
        graph = (graph0, graph1)
    except:
        pass

    try:
        layout0 = _read_h5_dict(filename, group="layout/layout_0")
        layout1 = _read_h5_dict(filename, group="layout/layout_1")
        layout = (layout0, layout1)
    except:
        pass

    try:
        germline = _read_h5_zip(
            filename, group="germline", key_group="keys", value_group="values"
        )
    except:

        pass

    try:
        with h5py.File(filename, "r") as hf:
            threshold = float(np.array(hf["threshold"]))
    except:
        threshold = None

    constructor = {}
    constructor["data"] = data
    if "metadata" in locals():
        constructor["metadata"] = metadata
    if "germline" in locals():
        constructor["germline"] = germline
    if "layout" in locals():
        constructor["layout"] = layout
    if "graph" in locals():
        constructor["graph"] = graph
    try:
        res = Dandelion(**constructor)
    except:
        res = Dandelion(**constructor, initialize=False)

    if "threshold" in locals():
        res.threshold = threshold
    else:
        pass
    return res


def read_airr(
    file: Path | str,
    prefix: str | None = None,
    suffix: str | None = None,
    sep: str = "_",
    remove_trailing_hyphen_number: bool = False,
) -> Dandelion:
    """
    Reads a standard single-cell AIRR rearrangement file.

    If you have non-single-cell data, use `.load_data` first to load the data and then pass it to `ddl.Dandelion`.
    That will tell you what columns are missing and you can fill it out accordingly e.g. make up a `cell_id`.

    Parameters
    ----------
    file : Path | str
        path to AIRR rearrangement .tsv file.
    prefix : str | None, optional
        Prefix to append to sequence_id and cell_id.
    suffix : str | None, optional
        Suffix to append to sequence_id and cell_id.
    sep : str, optional
        the separator to append suffix/prefix.
    remove_trailing_hyphen_number : bool, optional
        whether or not to remove the trailing hyphen number e.g. '-1' from the
        cell/contig barcodes.

    Returns
    -------
    Dandelion
        `Dandelion` object from AIRR file.
    """
    vdj = Dandelion(file)
    if suffix is not None:
        vdj.add_sequence_suffix(
            suffix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    elif prefix is not None:
        vdj.add_sequence_prefix(
            prefix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    return vdj


def read_bd_airr(
    file: Path | str,
    prefix: str | None = None,
    suffix: str | None = None,
    sep: str = "_",
    remove_trailing_hyphen_number: bool = False,
) -> Dandelion:
    """
    Read the TCR or BCR `_AIRR.tsv` produced from BD Rhapsody technology.

    Parameters
    ----------
    file : Path | str
        path to `_AIRR.tsv`
    prefix : str | None, optional
        Prefix to append to sequence_id and cell_id.
    suffix : str | None, optional
        Suffix to append to sequence_id and cell_id.
    sep : str, optional
        the separator to append suffix/prefix.
    remove_trailing_hyphen_number : bool, optional
        whether or not to remove the trailing hyphen number e.g. '-1' from the
        cell/contig barcodes.

    Returns
    -------
    Dandelion
        `Dandelion` object from BD AIRR file.
    """
    vdj = Dandelion(file)
    if suffix is not None:
        vdj.add_sequence_suffix(
            suffix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    elif prefix is not None:
        vdj.add_sequence_prefix(
            prefix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    return vdj


def read_parse_airr(
    file: Path | str,
    prefix: str | None = None,
    suffix: str | None = None,
    sep: str = "_",
    remove_trailing_hyphen_number: bool = False,
) -> Dandelion:
    """
    Read the TCR or BCR `_annotation_airr.tsv` produced from Parse Biosciences Evercode technology.

    This is not to be used for any airr rearrangement file, but specifically for the one produced by Parse Biosciences.
    For standard airr rearrangement files e.g. `all_contig_dandelion.tsv`, use `ddl.Dandelion` or `ddl.read_airr` directly.

    Parameters
    ----------
    file : Path | str
        path to `_annotation_airr.tsv`
    prefix : str | None, optional
        Prefix to append to sequence_id and cell_id.
    suffix : str | None, optional
        Suffix to append to sequence_id and cell_id.
    sep : str, optional
        the separator to append suffix/prefix.
    remove_trailing_hyphen_number : bool, optional
        whether or not to remove the trailing hyphen number e.g. '-1' from the
        cell/contig barcodes.

    Returns
    -------
    Dandelion
        `Dandelion` object from Parse AIRR file.
    """
    data = load_data(file)
    data.drop("cell_id", axis=1, inplace=True)  # it's the wrong cell_id
    data = data.rename(
        columns={
            "cell_barcode": "cell_id",
            "read_count": "consensus_count",
            "transcript_count": "umi_count",
            "cdr3": "junction",
            "cdr3_aa": "junction_aa",
        }
    )
    vdj = Dandelion(data)
    if suffix is not None:
        vdj.add_sequence_suffix(
            suffix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    elif prefix is not None:
        vdj.add_sequence_prefix(
            prefix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    return vdj


def read_10x_airr(
    file: Path | str,
    prefix: str | None = None,
    suffix: str | None = None,
    sep: str = "_",
    remove_trailing_hyphen_number: bool = False,
) -> Dandelion:
    """
    Read the `airr_rearrangement.tsv` produced from Cell Ranger directly and returns a `Dandelion` object.

    This is not to be used for any airr rearrangement file, but specifically for the one produced by 10x Genomics.
    For standard airr rearrangement files e.g. `all_contig_dandelion.tsv`, use `ddl.Dandelion` or `ddl.read_airr` directly.

    Parameters
    ----------
    file : Path | str
        path to `airr_rearrangement.tsv`
    prefix : str | None, optional
        Prefix to append to sequence_id and cell_id.
    suffix : str | None, optional
        Suffix to append to sequence_id and cell_id.
    sep : str, optional
        the separator to append suffix/prefix.
    remove_trailing_hyphen_number : bool, optional
        whether or not to remove the trailing hyphen number e.g. '-1' from the
        cell/contig barcodes.

    Returns
    -------
    Dandelion
        `Dandelion` object from 10x AIRR file.
    """
    dat = load_data(file)
    # get all the v,d,j,c calls
    if "locus" not in dat:
        tmp = [
            (v, d, j, c)
            for v, d, j, c in zip(
                dat["v_call"], dat["d_call"], dat["j_call"], dat["c_call"]
            )
        ]
        locus = []
        for t in tmp:
            if all("IGH" in x for x in t if pd.notnull(x)):
                locus.append("IGH")
            elif all("IGK" in x for x in t if pd.notnull(x)):
                locus.append("IGK")
            elif all("IGL" in x for x in t if pd.notnull(x)):
                locus.append("IGL")
            elif all("TRA" in x for x in t if pd.notnull(x)):
                locus.append("TRA")
            elif all("TRB" in x for x in t if pd.notnull(x)):
                locus.append("TRB")
            elif all("TRD" in x for x in t if pd.notnull(x)):
                locus.append("TRD")
            elif all("TRG" in x for x in t if pd.notnull(x)):
                locus.append("TRG")
            else:
                locus.append(np.nan)
        dat["locus"] = locus
    null_columns = [col for col in dat.columns if all_missing(dat[col])]
    if len(null_columns) > 0:
        dat.drop(null_columns, inplace=True, axis=1)
    vdj = Dandelion(dat)
    if suffix is not None:
        vdj.add_sequence_suffix(
            suffix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    elif prefix is not None:
        vdj.add_sequence_prefix(
            prefix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    return vdj


def read_10x_vdj(
    path: Path | str,
    filename_prefix: str | None = None,
    prefix: str | None = None,
    suffix: str | None = None,
    sep: str = "_",
    remove_malformed: bool = True,
    remove_trailing_hyphen_number: bool = False,
) -> Dandelion:
    """
    A parser to read .csv and .json files directly from folder containing 10x cellranger-outputs.

    This function parses the 10x output files into an AIRR compatible format.

    Minimum requirement is one of either {filename_prefix}_contig_annotations.csv or all_contig_annotations.json.

    If .fasta, .json files are found in the same folder, additional info will be appended to the final table.

    Parameters
    ----------
    path : Path | str
        path to folder containing `.csv` and/or `.json` files, or path to files directly.
    filename_prefix : str | None, optional
        prefix of file name preceding '_contig'. None defaults to 'all'.
    prefix : str | None, optional
        Prefix to append to sequence_id and cell_id.
    suffix : str | None, optional
        Suffix to append to sequence_id and cell_id.
    sep : str, optional
        the separator to append suffix/prefix.
    remove_malformed : bool, optional
        whether or not to remove malformed contigs.
    remove_trailing_hyphen_number : bool, optional
        whether or not to remove the trailing hyphen number e.g. '-1' from the
        cell/contig barcodes.

    Returns
    -------
    Dandelion
        `Dandelion` object holding the parsed data.

    Raises
    ------
    IOError
        if contig_annotations.csv and all_contig_annotations.json file(s) not found in the input folder.

    """
    filename_pre = (
        DEFAULT_PREFIX if filename_prefix is None else filename_prefix
    )

    if os.path.isdir(str(path)):
        files = os.listdir(path)
        filelist = []
        for fx in files:
            if re.search(filename_pre + "_contig", fx):
                if fx.endswith(".fasta") or fx.endswith(".csv"):
                    filelist.append(fx)
            if re.search(
                f"{filename_pre.replace('filtered', 'all')}_contig_annotations",
                fx,
            ):
                if fx.endswith(".json"):
                    filelist.append(fx)
        csv_idx = [i for i, j in enumerate(filelist) if j.endswith(".csv")]
        json_idx = [i for i, j in enumerate(filelist) if j.endswith(".json")]
        if len(csv_idx) == 1:
            file = str(path) + "/" + str(filelist[csv_idx[0]])
            logg.info("Reading {}".format(str(file)))
            raw = pd.read_csv(str(file))
            raw.set_index("contig_id", drop=False, inplace=True)
            fasta_file = str(file).split("_annotations.csv")[0] + ".fasta"
            json_file = re.sub(
                filename_pre + "_contig_annotations",
                f"{filename_pre.replace('filtered', 'all')}_contig_annotations",
                str(file).split(".csv")[0] + ".json",
            )
            if os.path.exists(json_file):
                logg.info(
                    "Found {} file. Extracting extra information.".format(
                        str(json_file)
                    )
                )
                out = parse_annotation(raw)
                with open(json_file) as f:
                    raw_json = json.load(f)
                out_json = parse_json(raw_json)
                out.update(out_json)
            elif os.path.exists(fasta_file):
                logg.info(
                    "Found {} file. Extracting extra information.".format(
                        str(fasta_file)
                    )
                )
                seqs = {}
                fh = open(fasta_file)
                for header, sequence in fasta_iterator(fh):
                    seqs[header] = sequence
                raw["sequence"] = pd.Series(seqs)
                out = parse_annotation(raw)
            else:
                out = parse_annotation(raw)
        elif len(csv_idx) < 1:
            if len(json_idx) == 1:
                json_file = str(path) + "/" + str(filelist[json_idx[0]])
                logg.info("Reading {}".format(json_file))
                if os.path.exists(json_file):
                    with open(json_file) as f:
                        raw = json.load(f)
                    out = parse_json(raw)
            else:
                raise OSError(
                    "{}_contig_annotations.csv and {}_contig_annotations.json file(s) not found in {} folder.".format(
                        str(filename_pre),
                        filename_pre.replace("filtered", "all"),
                        str(path),
                    )
                )
        elif len(csv_idx) > 1:
            raise OSError(
                "There are multiple input .csv files with the same filename prefix {} in {} folder.".format(
                    str(filename_pre), str(path)
                )
            )
    elif os.path.isfile(str(path)):
        file = path
        if str(file).endswith(".csv"):
            logg.info("Reading {}.".format(str(file)))
            raw = pd.read_csv(str(file))
            raw.set_index("contig_id", drop=False, inplace=True)
            fasta_file = str(file).split("_annotations.csv")[0] + ".fasta"
            json_file = re.sub(
                filename_pre + "_contig_annotations",
                f"{filename_pre.replace('filtered', 'all')}_contig_annotations",
                str(file).split(".csv")[0] + ".json",
            )
            if os.path.exists(json_file):
                logg.info(
                    "Found {} file. Extracting extra information.".format(
                        str(json_file)
                    )
                )
                out = parse_annotation(raw)
                with open(json_file) as f:
                    raw_json = json.load(f)
                out_json = parse_json(raw_json)
                out.update(out_json)
            elif os.path.exists(fasta_file):
                logg.info(
                    "Found {} file. Extracting extra information.".format(
                        str(fasta_file)
                    )
                )
                seqs = {}
                fh = open(fasta_file)
                for header, sequence in fasta_iterator(fh):
                    seqs[header] = sequence
                raw["sequence"] = pd.Series(seqs)
                out = parse_annotation(raw)
            else:
                out = parse_annotation(raw)
        elif str(file).endswith(".json"):
            if os.path.exists(file):
                logg.info("Reading {}".format(file))
                with open(file) as f:
                    raw = json.load(f)
                out = parse_json(raw)
            else:
                raise OSError("{} not found.".format(file))
    else:
        raise OSError("{} not found.".format(path))
    res = pd.DataFrame.from_dict(out, orient="index")
    # quick check if locus is malformed
    if remove_malformed:
        res = res[~res["locus"].str.contains("[|]")]
    vdj = Dandelion(res)
    if suffix is not None:
        vdj.add_sequence_suffix(
            suffix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    elif prefix is not None:
        vdj.add_sequence_prefix(
            prefix,
            sep=sep,
            remove_trailing_hyphen_number=remove_trailing_hyphen_number,
        )
    return vdj


def parse_json(data: list) -> defaultdict:
    """Parse json file."""
    main_dict1 = {
        "barcode": "cell_id",
        "contig_name": "sequence_id",
        "sequence": "sequence",
        "aa_sequence": "sequence_aa",
        "productive": "productive",
        "full_length": "complete_vdj",
        "frame": "vj_in_frame",
        "cdr3_seq": "junction",
        "cdr3": "junction_aa",
    }
    main_dict2 = {
        "read_count": "consensus_count",
        "umi_count": "umi_count",
        # "duplicate_count": "umi_count",
        "cdr3_start": "cdr3_start",
        "cdr3_stop": "cdr3_end",
    }
    main_dict3 = {
        "high_confidence": "high_confidence_10x",
        "filtered": "filtered_10x",
        "is_gex_cell": "is_cell_10x",
        "is_asm_cell": "is_asm_cell_10x",
    }
    info_dict = {
        "raw_clonotype_id": "clone_id",
        "raw_consensus_id": "raw_consensus_id_10x",
        "exact_subclonotype_id": "exact_subclonotype_id_10x",
    }
    region_type_dict = {
        "L-REGION+V-REGION": "v_call",
        "D-REGION": "d_call",
        "J-REGION": "j_call",
        "C-REGION": "c_call",
    }
    required_calls = ["v_call", "d_call", "j_call", "c_call"]
    region_keys = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4"]
    out = defaultdict(OrderedDict)
    for i in range(len(data)):
        if data[i]["contig_name"] is not None:
            key = data[i]["contig_name"]
        else:
            continue
        for k in main_dict1.keys():
            if k in data[i]:
                if data[i][k] is not None:
                    out[key].update({main_dict1[k]: data[i][k]})
                else:
                    out[key].update({main_dict1[k]: ""})
            else:
                out[key].update({main_dict1[k]: ""})
        if data[i]["annotations"] is not None:
            chains = []
            for j in range(len(data[i]["annotations"])):
                chains.append(data[i]["annotations"][j]["feature"]["chain"])
            out[key].update({"locus": "|".join(list(set(chains)))})
            for j in range(len(data[i]["annotations"])):
                rtype = data[i]["annotations"][j]["feature"]["region_type"]
                if rtype in region_type_dict.keys():
                    call = region_type_dict[rtype]
                else:
                    continue
                if (
                    data[i]["annotations"][j]["feature"]["gene_name"]
                    is not None
                ):
                    out[key].update(
                        {
                            call: data[i]["annotations"][j]["feature"][
                                "gene_name"
                            ]
                        }
                    )
                else:
                    out[key].update({call: ""})
        for rc in required_calls:
            if rc not in out[key]:
                out[key].update({rc: ""})
        for k in main_dict2.keys():
            if k in data[i]:
                if data[i][k] is not None:
                    out[key].update({main_dict2[k]: data[i][k]})
                else:
                    out[key].update({main_dict2[k]: np.nan})
            else:
                out[key].update({main_dict2[k]: np.nan})
        for rk in region_keys:
            if rk in data[i]:
                if data[i][rk] is not None:
                    for k in data[i][rk]:
                        if k == "start":
                            ka = rk + "_start"
                        elif k == "stop":
                            ka = rk + "_end"
                        elif k == "nt_seq":
                            ka = rk + ""
                        elif k == "aa_seq":
                            ka = rk + "_aa"
                        else:
                            continue
                        out[key].update({ka: data[i][rk][k]})
                else:
                    for k in region_keys:
                        out[key].update({k + "_start": np.nan})
                        out[key].update({k + "_end": np.nan})
                        out[key].update({k + "": ""})
                        out[key].update({k + "_aa": ""})
            else:
                for k in region_keys:
                    out[key].update({k + "_start": np.nan})
                    out[key].update({k + "_end": np.nan})
                    out[key].update({k + "": ""})
                    out[key].update({k + "_aa": ""})
        if data[i]["info"] is not None:
            for info in data[i]["info"]:
                if data[i]["info"][info] is not None:
                    out[key].update({info_dict[info]: data[i]["info"][info]})
                else:
                    out[key].update({info_dict[info]: ""})
        for k in main_dict3.keys():
            if k in data[i]:
                if data[i][k] is not None:
                    out[key].update({main_dict3[k]: data[i][k]})
                else:
                    out[key].update({main_dict3[k]: ""})
            else:
                out[key].update({main_dict3[k]: ""})
    return out


def parse_annotation(data: pd.DataFrame) -> defaultdict:
    """Parse annotation file."""
    out = defaultdict(OrderedDict)
    swap_dict = dict(zip(CELLRANGER, AIRR))
    for _, row in data.iterrows():
        contig = Contig(row, swap_dict).contig["sequence_id"]
        out[contig] = Contig(row, swap_dict).contig
        if out[contig]["locus"] in ["None", "none", None, np.nan, ""]:
            calls = []
            for call in ["v_call", "d_call", "j_call", "c_call"]:
                if out[contig][call] not in ["None", "none", None, np.nan, ""]:
                    calls.append(out[contig][call])
            out[contig]["locus"] = "|".join(list({str(c)[:3] for c in calls}))
        if out[contig]["locus"] == "None" or out[contig]["locus"] == "":
            out[contig]["locus"] = "|"
    return out


def change_file_location(
    data: list[Path | str],
    filename_prefix: list[str] | str | None = None,
) -> None:
    """
    Move file from tmp folder to dandelion folder.

    Only used for TCR data.

    Parameters
    ----------
    data : list[Path | str]
        list of data folders containing the .tsv files. if provided as a single string, it will first be converted to a
        list; this allows for the function to be run on single/multiple samples.
    filename_prefix : list[str] | str | None, optional
        list of prefixes of file names preceding '_contig'. None defaults to 'all'.
    """
    fileformat = "blast"
    data, filename_prefix = check_data(data, filename_prefix)

    informat_dict = {
        "changeo": "_igblast_db-pass.tsv",
        "blast": "_igblast_db-pass.tsv",
        "airr": "_igblast_gap.tsv",
    }
    filePath = None

    for i in range(0, len(data)):
        filePath = check_filepath(
            data[i],
            filename_prefix=filename_prefix[i],
            ends_with=informat_dict[fileformat],
            sub_dir="tmp",
        )
        if filePath is not None:
            tmp = check_travdv(filePath)
            _airrfile = str(filePath).replace("_db-pass.tsv", ".tsv")
            airr_output = load_data(_airrfile)
            cols_to_merge = [
                "junction_aa_length",
                "fwr1_aa",
                "fwr2_aa",
                "fwr3_aa",
                "fwr4_aa",
                "cdr1_aa",
                "cdr2_aa",
                "cdr3_aa",
                "sequence_alignment_aa",
                "v_sequence_alignment_aa",
                "d_sequence_alignment_aa",
                "j_sequence_alignment_aa",
            ]
            for x in cols_to_merge:
                tmp[x] = pd.Series(airr_output[x])

            write_airr(tmp, filePath)
            fp = Path(filePath)
            shutil.copyfile(fp, fp.parent.parent / fp.name)


def move_to_tmp(
    data: list[Path | str], filename_prefix: list[str] | str | None = None
) -> None:
    """Move file to tmp."""
    data, filename_prefix = check_data(data, filename_prefix)

    for i in range(0, len(data)):
        filePath1 = check_filepath(
            data[i],
            filename_prefix=filename_prefix[i],
            ends_with="_annotations.csv",
        )
        filePath2 = check_filepath(
            data[i], filename_prefix=filename_prefix[i], ends_with=".fasta"
        )
        for fp in [filePath1, filePath2]:
            fp = Path(fp)
            shutil.move(fp, fp.parent / "tmp" / fp.name)


def make_all(
    data: list[Path | str],
    filename_prefix: list[str] | str | None = None,
    loci: Literal["ig", "tr"] = "tr",
) -> None:
    """Construct db-all tsv file."""
    data, filename_prefix = check_data(data, filename_prefix)

    for i in range(0, len(data)):
        if loci == "tr":
            filePath1 = check_filepath(
                data[i],
                filename_prefix=filename_prefix[i],
                ends_with="_igblast_db-pass.tsv",
                sub_dir="tmp",
            )
        else:
            filePath1 = check_filepath(
                data[i],
                filename_prefix=filename_prefix[i],
                ends_with="_igblast_db-pass_genotyped.tsv",
                sub_dir="tmp",
            )
            if filePath1 is None:
                filePath1 = check_filepath(
                    data[i],
                    filename_prefix=filename_prefix[i],
                    ends_with="_igblast_db-pass.tsv",
                    sub_dir="tmp",
                )
                out_ex = "db-pass.tsv"
            else:
                out_ex = "db-pass_genotyped.tsv"
        filePath2 = check_filepath(
            data[i],
            filename_prefix=filename_prefix[i],
            ends_with="_igblast_db-fail.tsv",
            sub_dir="tmp",
        )
        if filePath1 is not None:
            df1 = pd.read_csv(filePath1, sep="\t")
            df1 = check_complete(df1)
            write_airr(df1, filePath1)
            if filePath2 is not None:
                df2 = pd.read_csv(filePath2, sep="\t")
                df2 = check_complete(df2)
                df = pd.concat([df1, df2])
                if loci == "tr":
                    write_airr(
                        df,
                        filePath1.parent
                        / (
                            filePath1.name.rsplit("db-pass.tsv")[0]
                            + "db-all.tsv"
                        ),
                    )
                else:
                    write_airr(
                        df,
                        filePath1.parent
                        / (filePath1.name.rsplit(out_ex)[0] + "db-all.tsv"),
                    )
                write_airr(df2, filePath2)
            else:
                if loci == "tr":
                    write_airr(
                        df1,
                        filePath1.parent
                        / (
                            filePath1.name.rsplit("db-pass.tsv")[0]
                            + "db-all.tsv"
                        ),
                    )
                else:
                    write_airr(
                        df1,
                        filePath1.parent
                        / (filePath1.name.rsplit(out_ex)[0] + "db-all.tsv"),
                    )


def rename_dandelion(
    data: list[Path | str],
    filename_prefix: list[str] | str | None = None,
    ends_with: str = "_igblast_db-pass_genotyped.tsv",
    sub_dir: str | None = None,
) -> None:
    """Rename final dandlion file."""
    data, filename_prefix = check_data(data, filename_prefix)

    for i in range(0, len(data)):
        filePath = check_filepath(
            data[i],
            filename_prefix=filename_prefix[i],
            ends_with=ends_with,
            sub_dir=sub_dir,
        )  # must be whatever's after contig
        if sub_dir is None:
            fp = filePath.parent / filePath.name.rsplit(ends_with)[0]
        else:
            fp = filePath.parent.parent / filePath.name.rsplit(ends_with)[0]
        shutil.move(filePath, Path(str(fp) + "_dandelion.tsv"))


def check_complete(df: pd.DataFrame) -> pd.DataFrame:
    """check if contig contains cdr3.

    Parameters
    ----------
    df : pd.DataFrame
        airr data frame.

    Returns
    -------
    pd.DataFrame
        completed airr data frame
    """
    if "complete_vdj" not in df:
        df["complete_vdj"] = ""
    for i in df.index:
        junc = df.loc[i, "junction"]
        if not present(junc):
            df.at[i, "productive"] = "F"
            df.at[i, "complete_vdj"] = "F"
    return df


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


def _read_h5_group(filename: Path | str, group: str) -> pd.DataFrame:
    """
    Read a specific group from an H5 file.

    Parameters
    ----------
    filename : Path | str
        The path to the H5 file.
    group : str
        The group to read from the H5 file.

    Returns
    -------
    pd.DataFrame
        The data from the specified group as a pandas dataframe.

    Raises
    ------
    KeyError
        If the specified group is not found in the H5 file.
    """
    try:
        with h5py.File(filename, "r") as hf:
            data_group = hf[group]
            structured_data_array = data_group[:]
            decoded = {}
            if structured_data_array.dtype.names is not None:
                for col in structured_data_array.dtype.names:
                    dtype = structured_data_array.dtype[col]
                    if dtype.char == "S":  # Check if dtype is byte strings
                        # Decode byte strings
                        decoded.update(
                            {
                                col: np.array(
                                    [
                                        (
                                            x.astype(str)
                                            if isinstance(x, bytes)
                                            else x
                                        )
                                        for x in structured_data_array[col]
                                    ]
                                )
                            }
                        )
                    else:
                        decoded.update({col: structured_data_array[col]})
                # Create a DataFrame from the structured array
                data = pd.DataFrame(decoded)
            else:
                data = np.array(
                    [
                        x.astype(str) if isinstance(x, bytes) else x
                        for x in structured_data_array
                    ]
                )
    except TypeError:
        try:
            data = pd.read_hdf(filename, group)
        except:
            raise KeyError(
                f"{str(filename)} does not contain attribute `{group}`"
            )
    return data


def _read_h5_csr_matrix(filename: Path | str, group: str) -> pd.DataFrame:
    """
    Read a group from an H5 file originally stored as a compressed sparse matrix.

    Parameters
    ----------
    filename : Path | str
        The path to the H5 file.
    group : str
        The group to read from the H5 file.

    Returns
    -------
    pd.DataFrame
        The data from the specified group as a pandas dataframe.
    """
    try:
        with h5py.File(filename, "r") as f:
            data = f[f"{group}/data"][:]
            indices = f[f"{group}/indices"][:]
            indptr = f[f"{group}/indptr"][:]
            shape = tuple(f[f"{group}/shape"][:])
            # Reconstruct CSR matrix
            loaded_matrix = csr_matrix((data, indices, indptr), shape=shape)
            df = pd.DataFrame(loaded_matrix.toarray())
            df_col = _read_h5_group(filename, f"{group}/column")
            df_index = _read_h5_group(filename, f"{group}/index")
            df.columns = df_col
            df.index = df_index
    except TypeError:
        try:
            data = pd.read_hdf(filename, group)
        except:
            raise KeyError(
                f"{str(filename)} does not contain attribute `{group}`"
            )
    return df


def _create_graph(
    adj: pd.DataFrame,
    adjust_adjacency: int | float = 0,
    fillna: int | float = 0,
) -> nx.Graph:
    """
    Create a networkx graph from the given adjacency matrix.

    Parameters
    ----------
    adj : pd.DataFrame
        The adjacency matrix to create the graph from.
    adjust_adjacency : int | float, optional
        The value to add to the graph by as a way to adjust the adjacency matrix. Defaults to 0.
    fillna : int | float, optional
        The value to fill NaN values with. Defaults to 0.

    Returns
    -------
    nx.Graph
        The created networkx graph.
    """
    if adjust_adjacency != 0:
        adj += adjust_adjacency
    adj = adj.fillna(fillna)
    g = nx.from_pandas_adjacency(adj)

    if adjust_adjacency != 0:
        for u, v, d in g.edges(data=True):
            d["weight"] -= adjust_adjacency

    return g


def _read_h5_dict(filename: Path | str, group: str) -> dict:
    """
    Read a dictionary from an H5 file.

    Parameters
    ----------
    filename : Path | str
        The path to the H5 file.
    group : str
        The group to read from the H5 file.

    Returns
    -------
    dict
        The dictionary data from the specified group.
    """
    out_dict = {}
    with h5py.File(filename, "r") as hf:
        for k in hf[group].keys():
            out_dict[k] = hf[group][k][:]

    return out_dict


def _read_h5_zip(
    filename: Path | str, group: str, key_group: str, value_group: str
) -> dict:
    """
    Read two groups from an H5 file and return them as a dictionary.

    Parameters
    ----------
    filename : Path | str
        The path to the H5 file.
    group : str
        The group to read from the H5 file.
    key_group : str
        The name of the group containing the keys.
    value_group : str
        The name of the group containing the values.

    Returns
    -------
    dict
        The dictionary data from the specified groups.
    """
    with h5py.File(filename, "r") as hf:
        try:
            keys = [
                key.decode("utf-8") for key in hf[f"{group}/{key_group}"][:]
            ]
            values = [
                value.decode("utf-8")
                for value in hf[f"{group}/{value_group}"][:]
            ]
            # Reconstruct the dictionary
            out_dict = dict(zip(keys, values))
        except:
            out_dict = {}
            for g in hf[group].attrs:
                out_dict.update({g: hf[group].attrs[g]})
    return out_dict
