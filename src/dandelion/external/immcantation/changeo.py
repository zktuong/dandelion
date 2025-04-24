from pathlib import Path
from scanpy import logging as logg
from subprocess import run
from typing import Literal

from dandelion.utilities._utilities import (
    set_germline_env,
    set_igblast_env,
    NO_DS,
)


def assigngenes_igblast(
    fasta: Path | str,
    igblast_db: Path | str | None = None,
    org: Literal["human", "mouse"] = "human",
    loci: Literal["ig", "tr"] = "ig",
    additional_args: list[str] = [],
):
    """
    Reannotate with IgBLASTn.

    Parameters
    ----------
    fasta : Path | str
        path to fasta file for reannotation.
    igblast_db : Path | str | None, optional
        path to igblast database.
    org : Literal["human", "mouse"], optional
        organism for germline sequences.
    loci : Literal["ig", "tr"], optional
        `ig` or `tr` mode for running igblastn.
    additional_args : list[str], optional
        Additional arguments to pass to `AssignGenes.py`.
    """
    env, igdb, fasta = set_igblast_env(igblast_db=igblast_db, input_file=fasta)
    outfolder = fasta.parent / "tmp"
    outfolder.mkdir(parents=True, exist_ok=True)

    informat_dict = {"blast": "_igblast.fmt7", "airr": "_igblast.tsv"}

    for fileformat in ["blast", "airr"]:
        outfile = fasta.stem + informat_dict[fileformat]
        cmd = [
            "AssignGenes.py",
            "igblast",
            "-s",
            str(fasta),
            "-b",
            str(igdb),
            "--organism",
            org,
            "--loci",
            loci,
            "--format",
            fileformat,
            "-o",
            str(outfolder / outfile),
        ]
        cmd += additional_args
        logg.info("Running command: %s\n" % (" ".join(cmd)))
        run(cmd, env=env)  # logs are printed to terminal


def makedb_igblast(
    fasta: Path | str,
    igblast_output: Path | str | None = None,
    germline: str | None = None,
    org: Literal["human", "mouse"] = "human",
    db: Literal["imgt", "ogrdb"] = "imgt",
    extended: bool = True,
    additional_args: list[str] = [],
    loci: Literal["ig", "tr"] = "ig",
):
    """
    Parse IgBLAST output to AIRR format.

    Parameters
    ----------
    fasta : Path | str
        path to fasta file used for reannotation.
    igblast_output : Path | str | None, optional
        path to igblast output file.
    germline : str | None, optional
        path to germline database.
    org : Literal["human", "mouse"], optional
        organism of germline sequences.
    db : Literal["imgt", "ogrdb"], optional
        `imgt` or `ogrdb` reference database for running igblastn.
    extended : bool, optional
        whether or not to parse extended 10x annotations.
    additional_args: list[str], optional
        Additional arguments to pass to `MakeDb.py`.
    """
    env, gml, fasta = set_germline_env(
        germline=germline,
        org=org,
        input_file=fasta,
        db=db,
    )
    if igblast_output is None:
        indir = fasta.parent / "tmp"
        infile = fasta.stem + "_igblast.fmt7"
        igbo = indir / infile
    else:
        igbo = Path(igblast_output)

    cellranger_annotation = fasta.parent / (fasta.stem + "_annotations.csv")

    if (org == "mouse") and (loci.lower() == "tr"):
        cmd = [
            "MakeDb_gentle.py",
            "igblast",
            "-i",
            str(igbo),
            "-s",
            str(fasta),
            "-r",
            str(gml),
            "--10x",
            str(cellranger_annotation),
        ]
    else:
        cmd = [
            "MakeDb.py",
            "igblast",
            "-i",
            str(igbo),
            "-s",
            str(fasta),
            "-r",
            str(gml),
            "--10x",
            str(cellranger_annotation),
        ]
    if extended:
        cmd = cmd + ["--extended"]
    for add_cmd in [[], ["--failed"]]:
        cmd = cmd + add_cmd + additional_args
        logg.info("Running command: %s\n" % (" ".join(cmd)))
        run(cmd, env=env)  # logs are printed to terminal


def parsedb_heavy(airr_file: Path | str):
    """
    Parse AIRR tsv file (heavy chain contigs only).

    Parameters
    ----------
    airr_file : Path | str
        path to AIRR tsv file.
    """
    outname = Path(airr_file).stem + "_heavy"
    cmd = [
        "ParseDb.py",
        "select",
        "-d",
        str(airr_file),
        "-f",
        "locus",
        "-u",
        "IGH",
        "--logic",
        "all",
        "--regex",
        "--outname",
        outname,
    ]

    logg.info("Running command: %s\n" % (" ".join(cmd)))
    run(cmd)  # logs are printed to terminal


def parsedb_light(airr_file: Path | str):
    """
    Parse AIRR tsv file (light chain contigs only).

    Parameters
    ----------
    airr_file : Path | str
        path to AIRR tsv file.
    """
    outname = Path(airr_file).stem + "_light"
    cmd = [
        "ParseDb.py",
        "select",
        "-d",
        str(airr_file),
        "-f",
        "locus",
        "-u",
        "IG[LK]",
        "--logic",
        "all",
        "--regex",
        "--outname",
        outname,
    ]

    logg.info("Running command: %s\n" % (" ".join(cmd)))
    run(cmd)  # logs are printed to terminal


def creategermlines(
    airr_file: Path | str,
    germline: list[str] | None = None,
    org: Literal["human", "mouse"] = "human",
    genotyped_fasta: str | None = None,
    mode: Literal["heavy", "light"] | None = None,
    db: Literal["imgt", "ogrdb"] = "imgt",
    strain: (
        Literal[
            "c57bl6",
            "balbc",
            "129S1_SvImJ",
            "AKR_J",
            "A_J",
            "BALB_c_ByJ",
            "BALB_c",
            "C3H_HeJ",
            "C57BL_6J",
            "C57BL_6",
            "CAST_EiJ",
            "CBA_J",
            "DBA_1J",
            "DBA_2J",
            "LEWES_EiJ",
            "MRL_MpJ",
            "MSM_MsJ",
            "NOD_ShiLtJ",
            "NOR_LtJ",
            "NZB_BlNJ",
            "PWD_PhJ",
            "SJL_J",
        ]
        | None
    ) = None,
    additional_args: list[str] = [],
):
    """
    Wrapper for CreateGermlines.py for reconstructing germline sequences.

    Parameters
    ----------
    airr_file : Path | str
        path to AIRR tsv file.
    germline : list[str] | None, optional
        location to germline fasta files as a list.
    org : Literal["human", "mouse"], optional
        organism for germline sequences.
    genotyped_fasta : str | None, optional
        location to V genotyped fasta file.
    mode : Literal["heavy", "light"] | None, optional
        whether to run on heavy or light mode. If left as None, heavy and
        light will be run together.
    db : Literal["imgt", "ogrdb"], optional
        `imgt` or `ogrdb` reference database.
    strain : Literal["c57bl6", "balbc", "129S1_SvImJ", "AKR_J", "A_J", "BALB_c_ByJ", "BALB_c", "C3H_HeJ", "C57BL_6J", "C57BL_6", "CAST_EiJ", "CBA_J", "DBA_1J", "DBA_2J", "LEWES_EiJ", "MRL_MpJ", "MSM_MsJ", "NOD_ShiLtJ", "NOR_LtJ", "NZB_BlNJ", "PWD_PhJ", "SJL_J"] | None, optional
        strain of mouse to use for germline sequences. Only for `db="ogrdb"`. Note that only "c57bl6", "balbc", "CAST_EiJ", "LEWES_EiJ", "MSM_MsJ", "NOD_ShiLt_J" and "PWD_PhJ" contains both heavy chain and light chain germline sequences as a set.
        The rest will not allow igblastn and MakeDB.py to generate a successful airr table (check the failed file). "c57bl6" and "balbc" are merged databases of "C57BL_6" with "C57BL_6J" and "BALB_c" with "BALB_c_ByJ" respectively. None defaults to all combined.
    additional_args : list[str], optional
        Additional arguments to pass to `CreateGermlines.py`.
    """
    env, gml, airr_file = set_germline_env(
        germline=germline,
        org=org,
        input_file=airr_file,
        db=db,
    )
    _strain = "_" + strain if strain is not None else ""
    if germline is None:
        if mode == "heavy":
            if genotyped_fasta is None:
                gml_ref = [
                    str(gml / ("imgt_" + org + _strain + "_IGHV.fasta")),
                ]
            else:
                gml_ref = [
                    str(genotyped_fasta),
                ]
            gml_ref += [
                str(gml / ("imgt_" + org + _strain + "_IGHJ.fasta")),
            ]
            _strainD = "" if strain not in NO_DS else _strain
            gml_ref += [
                str(gml / ("imgt_" + org + _strainD + "_IGHD.fasta")),
            ]
        elif mode == "light":
            gml_ref = [
                str(gml / ("imgt_" + org + _strain + "_IGKV.fasta")),
                str(gml / ("imgt_" + org + _strain + "_IGKJ.fasta")),
                str(gml / ("imgt_" + org + _strain + "_IGLV.fasta")),
                str(gml / ("imgt_" + org + _strain + "_IGLJ.fasta")),
            ]
        elif mode is None:
            if genotyped_fasta is None:
                gml_ref = [
                    str(gml / ("imgt_" + org + _strain + "_IGHV.fasta")),
                ]
            else:
                gml_ref = [
                    str(genotyped_fasta),
                ]
            gml_ref += [
                str(gml / ("imgt_" + org + _strain + "_IGHJ.fasta")),
            ]
            _strainD = "" if strain not in NO_DS else _strain
            gml_ref += [
                str(gml / ("imgt_" + org + _strainD + "_IGHD.fasta")),
            ]
            gml_ref += [
                str(gml / ("imgt_" + org + _strain + "_IGKV.fasta")),
                str(gml / ("imgt_" + org + _strain + "_IGKJ.fasta")),
                str(gml / ("imgt_" + org + _strain + "_IGLV.fasta")),
                str(gml / ("imgt_" + org + _strain + "_IGLJ.fasta")),
            ]
    else:
        if not isinstance(germline, list):
            germline = [str(germline)]
        gml_ref = germline
    cmd = [
        "CreateGermlines.py",
        "-d",
        str(airr_file),
        "-r",
    ]
    cmd = cmd + gml_ref + additional_args

    logg.info("Running command: %s\n" % (" ".join(cmd)))
    run(cmd, env=env)  # logs are printed to terminal
