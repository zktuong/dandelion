#!/usr/bin/env python
import h5py
import os
import re
import warnings

import numpy as np
import pandas as pd

from airr import RearrangementSchema
from collections import defaultdict
from pathlib import Path
from subprocess import run
from typing import (
    TypeVar,
    Literal,
    Callable,
)

# help silence the dtype warning?
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

F = TypeVar("F", bound=Callable)  # Define a TypeVar for any callable type
MuData = TypeVar("mudata._core.mudata.MuData")
PResults = TypeVar("palantir.presults.PResults")

TRUES = ["T", "True", "true", "TRUE", True]
FALSES = ["F", "False", "false", "FALSE", False]
HEAVYLONG = ["IGH", "TRB", "TRD"]
LIGHTSHORT = ["IGK", "IGL", "TRA", "TRG"]
VCALL = "v_call"
JCALL = "j_call"
VCALLG = "v_call_genotyped"
JCALLG = "j_call_genotyped"
STRIPALLELENUM = "[*][0-9][0-9]"
NO_DS = [
    "129S1_SvImJ",
    "AKR_J",
    "A_J",
    "C3H_HeJ",
    "C57BL_6J",
    "BALB_c_ByJ",
    "CBA_J",
    "DBA_1J",
    "DBA_2J",
    "MRL_MpJ",
    "NOR_LtJ",
    "NZB_BlNJ",
    "SJL_J",
]
EMPTIES = [
    None,
    np.nan,
    pd.NA,
    "nan",
    "NaN",
    "",
]
DEFAULT_PREFIX = "all"

# for compatibility with python>=3.10
try:
    from collections.abc import Iterable
except ImportError:
    from collections.abc import Iterable


class Tree(defaultdict):
    """Create a recursive defaultdict."""

    def __init__(self, value=None) -> None:
        super().__init__(Tree)
        self.value = value


def dict_from_table(meta: pd.DataFrame, columns: tuple[str, str]) -> dict:
    """
    Generate a dictionary from a dataframe.

    Parameters
    ----------
    meta : pd.DataFrame
        pandas data frame or file path
    columns : tuple[str, str]
        column names in data frame

    Returns
    -------
    dict
        dictionary
    """
    if (isinstance(meta, pd.DataFrame)) & (columns is not None):
        meta_ = meta
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))
    elif (os.path.isfile(str(meta))) & (columns is not None):
        meta_ = pd.read_csv(meta, sep="\t", dtype="object")
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))

    sample_dict = clean_nan_dict(sample_dict)
    return sample_dict


def clean_nan_dict(d: dict) -> dict:
    """
    Remove nan from dictionary.

    Parameters
    ----------
    d : dict
        dictionary

    Returns
    -------
    dict
        dictionary with no NAs.
    """
    return {k: v for k, v in d.items() if v is not np.nan}


def flatten(l: list) -> list:
    """
    Flatten a list-in-list-in-list.

    Parameters
    ----------
    l : list
        a list-in-list list

    Yields
    ------
    list
        a flattened list.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def makeblastdb(ref: Path | str) -> None:
    """
    Run makeblastdb on constant region fasta file.

    Wrapper for makeblastdb.

    Parameters
    ----------
    ref : Path | str
        constant region fasta file.
    """
    cmd = ["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str(ref)]
    run(cmd)


def bh(pvalues: np.array) -> np.array:
    """
    Compute the Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    pvalues : np.array
        array of p-values to correct

    Returns
    -------
    np.array
        np.array of corrected p-values
    """
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n / rank) * pvalue)
    for i in range(0, int(n) - 1):
        if new_values[i] < new_values[i + 1]:
            new_values[i + 1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def is_categorical(array_like: pd.Series) -> bool:
    """Check if categorical."""
    return array_like.dtype.name == "category"


def type_check(dataframe: pd.DataFrame, key: str) -> bool:
    """Check dtype."""
    return (
        dataframe[key].dtype == str
        or dataframe[key].dtype == object
        or is_categorical(dataframe[key])
        or dataframe[key].dtype == bool
    )


def isGZIP(filename: str) -> bool:
    """Check if is gzipped file."""
    if filename.split(".")[-1] == "gz":
        return True
    return False


def isBZIP(filename: str) -> bool:
    """Check if is bzipped file."""
    if filename.split(".")[-1] == "pbz2":
        return True
    return False


def check_filepath(
    file_or_folder_path: Path | str,
    filename_prefix: str | None = None,
    ends_with: str | None = None,
    sub_dir: str | None = None,
    within_dandelion: bool = True,
) -> Path | None:
    """
    Checks whether file path exists.

    Parameters
    ----------
    file_or_folder_path : Path | str
        either a string or Path object pointing to a file or folder.
    filename_prefix : str | None, optional
        the prefix of the filename.
    ends_with : str | None, optional
        the suffix of the filename. Can be flexible i.e. not just the extension.
    sub_dir : str | None, optional
        the subdirectory to look for the file if specified
    within_dandelion : bool, optional
        whether to look for the file within a 'dandelion' sub folder.

    Returns
    -------
    Path | None
        Path object if file is found, else None.
    """
    filename_pre = (
        DEFAULT_PREFIX if filename_prefix is None else filename_prefix
    )

    ends_with = "" if ends_with is None else ends_with
    input_path = (
        Path(str(file_or_folder_path)).expanduser()
        if str(file_or_folder_path)[0] == "~"
        else Path(str(file_or_folder_path))
    )
    if input_path.is_file() and str(input_path).endswith(ends_with):
        return input_path
    elif input_path.is_dir():
        if within_dandelion:
            for child in input_path.iterdir():
                if child.name[0] != ".":
                    if child.is_dir() and child.name == "dandelion":
                        out_dir = child
                        if sub_dir is not None:
                            out_dir = out_dir / sub_dir
                        for file in out_dir.iterdir():
                            if file.name[0] != ".":
                                if file.is_file() and str(file).endswith(
                                    ends_with
                                ):
                                    if file.name.startswith(
                                        filename_pre + "_contig"
                                    ):
                                        return file
        else:
            if sub_dir is not None:
                input_path = input_path / sub_dir
            for file in input_path.iterdir():
                if file.name[0] != ".":
                    if file.is_file() and str(file).endswith(ends_with):
                        if file.name.startswith(filename_pre + "_contig"):
                            return file
    else:
        return None


def cmp(a, b):
    """Python2.x cmp function."""
    return (a > b) - (a < b)


def cmp_str_emptylast(s1, s2):
    """Help sort empty string to last."""
    if not s1 or not s2:
        return bool(s2) - bool(s1)

    return cmp(s1, s2)


def cmp_to_key(mycmp):
    """Convert a cmp= function into a key= function."""

    class K:
        """Key class"""

        def __init__(self, obj, *args) -> None:
            self.obj = obj

        def __lt__(self, other) -> bool:
            """Less than."""
            return mycmp(self.obj, other.obj) < 0

        def __gt__(self, other) -> bool:
            """Greater than."""
            return mycmp(self.obj, other.obj) > 0

        def __eq__(self, other) -> bool:
            """Equal."""
            return mycmp(self.obj, other.obj) == 0

        def __le__(self, other) -> bool:
            """Less than or equal."""
            return mycmp(self.obj, other.obj) <= 0

        def __ge__(self, other) -> bool:
            """Greater than or equal."""
            return mycmp(self.obj, other.obj) >= 0

        def __ne__(self, other) -> bool:
            """Not equal."""
            return mycmp(self.obj, other.obj) != 0

    return K


def not_same_call(a: str, b: str, pattern: str) -> bool:
    """Utility function to check if a == b in terms of pattern."""
    return (re.search(pattern, a) and not re.search(pattern, b)) or (
        re.search(pattern, b) and not re.search(pattern, a)
    )


def same_call(a: str, b: str, c: str, pattern: str) -> bool:
    """Utility function to check if a == b == c in terms of pattern."""
    queries = [a, b, c]
    queries = [q for q in queries if pd.notnull(q)]
    return all([re.search(pattern, x) for x in queries])


def present(x: str | None) -> bool:
    """Utility function to check if x is not null or blank."""
    return pd.notnull(x) and x not in [
        "",
        "None",
        "none",
        "NA",
        "na",
        "NaN",
        "nan",
    ]


def check_missing(x: str | None) -> bool:
    """Utility function to check if x is null or blank."""
    return pd.isnull(x) or x == ""


def all_missing(x: str | None) -> bool:
    """Utility function to check if all x is not null or blank."""
    return all(pd.isnull(x)) or all(x == "")


def all_missing2(x: str | None) -> bool:
    """Utility function to check if all x is not null or blank or the word None."""
    return all(pd.isnull(x)) or all(x == "") or all(x == "None")


def return_mix_dtype(data: pd.DataFrame) -> list:
    """Utility function to return mixed dtypes columns."""
    check = [
        c
        for c in data.columns
        if pd.api.types.infer_dtype(data[c]).startswith("mixed")
    ]
    return check


def get_numpy_dtype(series: pd.Series) -> str:
    """
    Map a Pandas dtype to an appropriate NumPy dtype.

    Parameters
    ----------
    series : pd.Series
        The Pandas Series.

    Returns
    -------
    str
        A string representing the NumPy dtype corresponding to the Pandas dtype.

    Raises
    ------
    TypeError
        If the Pandas dtype is unsupported.
    """
    if pd.api.types.is_integer_dtype(series):
        return "i4"  # 32-bit integer
    elif pd.api.types.is_float_dtype(series):
        return "f8"  # 64-bit float
    elif pd.api.types.is_bool_dtype(series):
        return "i1"  # 8-bit integer for booleans (True/False)
    elif pd.api.types.is_string_dtype(series) or pd.api.types.is_object_dtype(
        series
    ):
        # Handle object or string columns; dynamically calculate the max string length
        max_length = series.astype(str).map(len).max()
        return "S{}".format(max(1, max_length))  # String with max length
    else:
        raise TypeError(f"Unsupported data type: {series.name}")


def sanitize_data_for_saving(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """
    Quick sanitize dtypes for saving.

    Parameters
    ----------
    data : pd.DataFrame
        input dataframe.

    Returns
    -------
    tuple[pd.DataFrame, dict[str, str]]
        dataframe and corresponding numpy structured dtype.
    """
    tmp = data.copy()
    dtype_dict = {}
    for col in tmp:
        if col in RearrangementSchema.properties:
            dtype = RearrangementSchema.properties[col]["type"]
            tmp[col] = sanitize_column(tmp[col], dtype)
        else:
            tmp[col] = try_numeric_conversion(tmp[col])
        dtype_dict[col] = get_numpy_dtype(tmp[col])
    dtypes = [(key, record) for key, record in dtype_dict.items()]
    return tmp, dtypes


def sanitize_column(series: pd.Series, dtype: str) -> pd.Series:
    """
    Sanitize a column based on the specified dtype.

    Parameters
    ----------
    series : pd.Series
        The column to be sanitized.
    dtype : str
        The expected data type of the column (`string`, `boolean`, `integer`, or `number`).

    Returns
    -------
    pd.Series
        The sanitized column with replaced values and appropriate data type.
    """
    if dtype in ["string", "boolean"]:
        series = series.apply(lambda x: "" if pd.isna(x) else x)
        series = series.replace([None, np.nan, "nan", "na", "NaN", ""], "")
        return series.astype(str)
    elif dtype in ["integer", "number"]:
        series = series.apply(lambda x: np.nan if pd.isna(x) else x)
        return series.replace([None, np.nan, "nan", "na", "NaN", ""], np.nan)
    return series


def try_numeric_conversion(series: pd.Series) -> pd.Series:
    """
    Attempt to convert a column to numeric, or fallback to treating it as a string.

    Parameters
    ----------
    series : pd.Series
        The column to be converted.

    Returns
    -------
    pd.Series
        The column converted to numeric if possible, or sanitized as a string if not.
    """
    if series.dtype.name == "category":
        series = sanitize_column(series, "string")

    if series.apply(lambda x: isinstance(x, str) and "|" in x).any():
        return sanitize_column(series, "string")
    try:
        return pd.to_numeric(series)
    except:
        return sanitize_column(series, "string")


def sanitize_data(data: pd.DataFrame, ignore: str = "clone_id") -> None:
    """Quick sanitize dtypes."""
    data = data.astype("object")
    data = data.infer_objects()
    for d in data:
        if d in RearrangementSchema.properties:
            if RearrangementSchema.properties[d]["type"] in [
                "string",
                "boolean",
                "integer",
            ]:
                data[d] = data[d].replace(
                    EMPTIES,
                    "",
                )
                if RearrangementSchema.properties[d]["type"] == "integer":
                    data[d] = [
                        int(x) if present(x) else ""
                        for x in pd.to_numeric(data[d])
                    ]
            else:
                data[d] = data[d].replace(
                    EMPTIES,
                    np.nan,
                )
        else:
            if d != ignore:
                try:
                    data[d] = pd.to_numeric(data[d])
                except:
                    data[d] = data[d].replace(
                        to_replace=EMPTIES,
                        value="",
                    )
        if re.search("mu_freq", d):
            data[d] = [
                float(x) if present(x) else np.nan
                for x in pd.to_numeric(data[d])
            ]
        if re.search("mu_count", d):
            data[d] = [
                int(x) if present(x) else "" for x in pd.to_numeric(data[d])
            ]
    try:
        data = check_travdv(data)
    except:
        pass

    if (
        pd.Series(["cell_id", "umi_count", "productive"])
        .isin(data.columns)
        .all()
    ):  # sort so that the productive contig with the largest umi is first
        data.sort_values(
            by=["cell_id", "productive", "umi_count"],
            inplace=True,
            ascending=[True, False, False],
        )

    # check if airr-standards is happy
    validate_airr(data)
    return data


def sanitize_blastn(data: pd.DataFrame) -> None:
    """Quick sanitize dtypes."""
    data = data.astype("object")
    data = data.infer_objects()
    for d in data:
        if d in RearrangementSchema.properties:
            if RearrangementSchema.properties[d]["type"] in [
                "string",
                "boolean",
                "integer",
            ]:
                data[d] = data[d].replace(
                    EMPTIES,
                    "",
                )
                if RearrangementSchema.properties[d]["type"] == "integer":
                    data[d] = [
                        int(x) if present(x) else ""
                        for x in pd.to_numeric(data[d])
                    ]
            else:
                data[d] = data[d].replace(
                    EMPTIES,
                    np.nan,
                )
        else:
            try:
                data[d] = pd.to_numeric(data[d])
            except:
                data[d] = data[d].replace(
                    to_replace=EMPTIES,
                    value="",
                )
    return data


def validate_airr(data: pd.DataFrame) -> None:
    """Validate dtypes in airr table."""
    tmp = data.copy()
    int_columns = []
    for d in tmp:
        try:
            tmp[d].replace(np.nan, pd.NA).astype("Int64")
            int_columns.append(d)
        except:
            pass
    for _, row in tmp.iterrows():
        contig = Contig(row).contig
        for required in [
            "sequence",
            "rev_comp",
            "sequence_alignment",
            "germline_alignment",
            "v_cigar",
            "d_cigar",
            "j_cigar",
        ]:
            if required not in contig:
                contig.update({required: ""})
    RearrangementSchema.validate_header(contig.keys())
    RearrangementSchema.validate_row(contig)


def check_travdv(data: pd.DataFrame) -> pd.DataFrame:
    """Check if locus is TRA/D."""
    data = load_data(data)
    contig = [x for x in data["sequence_id"]]
    v = [x for x in data["v_call"]]
    d = [x for x in data["d_call"]]
    j = [x for x in data["j_call"]]
    c = [x for x in data["c_call"]]
    l = [x for x in data["locus"]]
    v_dict = dict(zip(contig, v))
    d_dict = dict(zip(contig, d))
    j_dict = dict(zip(contig, j))
    c_dict = dict(zip(contig, c))
    l_dict = dict(zip(contig, l))
    for co in contig:
        if re.search("TRAV.*/DV", v_dict[co]):
            if same_call(j_dict[co], c_dict[co], d_dict[co], "TRA"):
                if not re.search("TRA", l_dict[co]):
                    l_dict[co] = "TRA"
            elif same_call(j_dict[co], c_dict[co], d_dict[co], "TRD"):
                if not re.search("TRD", l_dict[co]):
                    l_dict[co] = "TRD"
    data["locus"] = pd.Series(l_dict)
    return data


def load_data(obj: pd.DataFrame | Path | str | None) -> pd.DataFrame:
    """
    Read in or copy dataframe object and set sequence_id as index without dropping.

    Parameters
    ----------
    obj : pd.DataFrame | Path | str | None
        file path to .tsv file or pandas DataFrame object.

    Returns
    -------
    pd.DataFrame

    Raises
    ------
    FileNotFoundError
        if input is not found.
    KeyError
        if `sequence_id` not found in input.
    """
    if obj is not None:
        if os.path.isfile(str(obj)):
            try:
                obj_ = pd.read_csv(obj, sep="\t")
            except FileNotFoundError as e:
                print(e)
        elif isinstance(obj, pd.DataFrame):
            obj_ = obj.copy()
        else:
            raise FileNotFoundError(
                "Either input is not of <class 'pandas.core.frame.DataFrame'> or file does not exist."
            )

        if "sequence_id" in obj_.columns:
            # assert that sequence_id is string
            obj_["sequence_id"] = obj_["sequence_id"].astype(str)
            obj_.set_index("sequence_id", drop=False, inplace=True)
            if "cell_id" not in obj_.columns:
                obj_["cell_id"] = [
                    c.split("_contig")[0] for c in obj_["sequence_id"]
                ]
            # assert that cell_id is string
            obj_["cell_id"] = obj_["cell_id"].astype(str)
        else:
            raise KeyError("'sequence_id' not found in columns of input")

        if "duplicate_count" in obj_.columns:
            if "umi_count" not in obj_.columns:
                obj_.rename(
                    columns={"duplicate_count": "umi_count"}, inplace=True
                )

        return obj_


class ContigDict(dict):
    """Class Object to extract the contigs as a dictionary."""

    def __setitem__(self, key: str, value: str) -> None:
        """Standard __setitem__."""
        super().__setitem__(key, value)

    def __hash__(self) -> int:
        """Make it hashable."""
        return hash(tuple(self))


class Contig:
    """Class Object to hold contig."""

    def __init__(self, contig: dict, mapper: dict | None = None) -> None:
        if mapper is not None:
            mapper.update({k: k for k in contig.keys() if k not in mapper})
            self._contig = ContigDict(
                {mapper[key]: vals for (key, vals) in contig.items()}
            )
        else:
            self._contig = ContigDict(contig)
        for key, value in self._contig.items():
            if isinstance(value, float) and np.isnan(value):
                self._contig[key] = ""

    @property
    def contig(self) -> "ContigDict":
        """Contig slot."""
        return self._contig


def mask_dj(
    data: list[Path | str],
    filename_prefix: str,
    d_evalue_threshold: float,
    j_evalue_threshold: float,
) -> None:
    """Mask d/j assignment."""
    for i in range(0, len(data)):
        filePath = check_filepath(
            data[i],
            filename_prefix=filename_prefix[i],
            ends_with="_igblast_db-pass.tsv",
        )
        if filePath is not None:
            dat = load_data(filePath)
            if "d_support_blastn" in dat:
                dat["d_call"] = [
                    "" if s > d_evalue_threshold else c
                    for c, s in zip(dat["d_call"], dat["d_support_blastn"])
                ]
            if "j_support_blastn" in dat:
                dat["j_call"] = [
                    "" if s > j_evalue_threshold else c
                    for c, s in zip(dat["j_call"], dat["j_support_blastn"])
                ]

            write_airr(dat, filePath)


def write_airr(data: pd.DataFrame, save: Path | str) -> None:
    """Save as airr formatted file."""
    data = sanitize_data(data)
    data.to_csv(save, sep="\t", index=False)


def write_blastn(data: pd.DataFrame, save: Path | str) -> None:
    """Write blast output."""
    data = sanitize_blastn(data)
    data.to_csv(save, sep="\t", index=False)


## from skbio==0.5.6
def _validate_counts_vector(
    counts: np.array, suppress_cast: bool = False
) -> np.array:
    """Validate and convert input to an acceptable counts vector type.
    Note: may not always return a copy of `counts`!
    """
    counts = np.asarray(counts)
    if not np.all(np.isreal(counts)):
        raise ValueError("Counts vector must contain real-valued entries.")
    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")
    elif (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts


def deprecated(
    details: str, deprecated_in: str, removed_in: str
) -> Callable[[F], F]:
    """Decorator to mark a function as deprecated"""

    def deprecated_decorator(func: F) -> F:
        """Deprecate dectorator"""

        def deprecated_func(*args, **kwargs):
            """Deprecate function"""
            warnings.warn(
                "{} is a deprecated in {} and will be removed in {}."
                " {}".format(func.__name__, deprecated_in, removed_in, details),
                category=DeprecationWarning,
                stacklevel=2,
            )
            return func(*args, **kwargs)

        return deprecated_func

    return deprecated_decorator


def format_call(
    metadata: pd.DataFrame,
    call: str,
    suffix_vdj: str = "_VDJ",
    suffix_vj: str = "_VJ",
) -> list:
    """Extract v/d/j/c call values from data."""
    call_dict = {
        "Multi": "Multi",
        "None": "None",
        "": "None",
        "unassigned": "None",
    }
    if suffix_vj is not None:
        call_1 = {
            x[0]: x[1] if present(x[1]) else "None"
            for x, y in zip(
                metadata[call + suffix_vdj].items(),
                list(metadata[call + suffix_vj]),
            )
        }
        call_2 = {
            x[0]: x[1] if present(x[1]) else "None"
            for x, y in zip(
                metadata[call + suffix_vj].items(),
                list(metadata[call + suffix_vdj]),
            )
        }
        call_2 = {x: y if "|" not in y else "Multi" for x, y in call_2.items()}
        call_4 = {
            x: "Single" if y not in call_dict else call_dict[y]
            for x, y in call_2.items()
        }
    else:
        call_1 = {
            x: y if present(y) else "None"
            for x, y in metadata[call + suffix_vdj].items()
        }
        call_2 = {x: "None" for x in call_1.keys()}
        call_4 = call_3 = call_2
    call_1 = {x: y if "|" not in y else "Multi" for x, y in call_1.items()}
    call_3 = {
        x: "Single" if y not in call_dict else call_dict[y]
        for x, y in call_1.items()
    }
    return (
        list(call_1.values()),
        list(call_2.values()),
        list(call_3.values()),
        list(call_4.values()),
    )


def format_locus(
    metadata: pd.DataFrame,
    suffix_vdj: str = "_VDJ",
    suffix_vj: str = "_VJ",
    productive_only: bool = True,
) -> pd.Series:
    """Extract locus call value from data."""
    locus_1 = dict(metadata["locus" + suffix_vdj])
    locus_2 = dict(metadata["locus" + suffix_vj])
    constant_1 = dict(metadata["isotype_status"])
    prod_1 = dict(metadata["productive" + suffix_vdj])
    prod_2 = dict(metadata["productive" + suffix_vj])
    locus_dict = {}
    for i in metadata.index:
        if productive_only:
            loc1 = {
                e: l
                for e, l in enumerate(
                    [
                        ll
                        for ll, p in zip(
                            locus_1[i].split("|"), prod_1[i].split("|")
                        )
                        if p in TRUES
                    ]
                )
            }
            loc2 = {
                e: l
                for e, l in enumerate(
                    [
                        ll
                        for ll, p in zip(
                            locus_2[i].split("|"), prod_2[i].split("|")
                        )
                        if p in TRUES
                    ]
                )
            }
        else:
            loc1 = {
                e: l for e, l in enumerate([ll for ll in locus_1[i].split("|")])
            }
            loc2 = {
                e: l for e, l in enumerate([ll for ll in locus_2[i].split("|")])
            }
        loc1x, loc2x = [], []
        if not all([px == "None" for px in loc1.values()]):
            loc1xx = list(loc1.values())
            loc1x = [ij[:2] for ij in loc1.values()]

        if not all([px == "None" for px in loc2.values()]):
            loc2xx = list(loc2.values())
            loc2x = [ij[:2] for ij in loc2.values()]

        if len(loc1x) > 0:
            if len(list(set(loc1x))) > 1:
                tmp1 = "ambiguous"
                if len(loc2x) > 0:
                    if len(list(set(loc2x))) > 1:
                        tmp2 = "ambiguous"
                    else:
                        if len(loc2x) > 1:
                            if (all(x in ["TRA", "TRG"] for x in loc2xx)) and (
                                len(list(set(loc2xx))) == 2
                            ):
                                tmp2 = "Extra VJ-exception"
                            else:
                                tmp2 = "Extra VJ"
                        else:
                            tmp2 = loc2xx[0]
                else:
                    tmp2 = "None"
            else:
                if len(loc1x) > 1:
                    if constant_1[i] == "IgM/IgD":
                        tmp1 = "IgM/IgD"
                    elif (all(x in ["TRB", "TRD"] for x in loc1xx)) and (
                        len(list(set(loc1xx))) == 2
                    ):
                        tmp1 = "Extra VDJ-exception"
                    else:
                        tmp1 = "Extra VDJ"
                else:
                    tmp1 = loc1xx[0]

                if len(loc2x) > 0:
                    if len(list(set(loc2x))) > 1:
                        tmp2 = "ambiguous"
                    else:
                        if len(loc2x) > 1:
                            if (all(x in ["TRA", "TRG"] for x in loc2xx)) and (
                                len(list(set(loc2xx))) == 2
                            ):
                                tmp2 = "Extra VJ-exception"
                            else:
                                tmp2 = "Extra VJ"
                        else:
                            tmp2 = loc2xx[0]
                else:
                    tmp2 = "None"

                if (
                    tmp1 not in ["None", "Extra VDJ", "Extra VDJ-exception"]
                ) and (tmp2 not in ["None", "Extra VJ", "Extra VJ-exception"]):
                    if list(set(loc1x)) != list(set(loc2x)):
                        tmp1 = "ambiguous"
                        tmp2 = "ambiguous"
        else:
            tmp1 = "None"
            if len(loc2x) > 0:
                if len(list(set(loc2x))) > 1:
                    tmp2 = "ambiguous"
                else:
                    if len(loc2x) > 1:
                        if (all(x in ["TRA", "TRG"] for x in loc2xx)) and (
                            len(list(set(loc2xx))) == 2
                        ):
                            tmp2 = "Extra VJ-exception"
                        else:
                            tmp2 = "Extra VJ"
                    else:
                        tmp2 = loc2xx[0]
            else:
                tmp2 = "None"
        if any(tmp == "ambiguous" for tmp in [tmp1, tmp2]):
            locus_dict.update({i: "ambiguous"})
        else:
            locus_dict.update({i: tmp1 + " + " + tmp2})

        if any(tmp == "None" for tmp in [tmp1, tmp2]):
            if tmp1 == "None":
                locus_dict.update({i: "Orphan " + tmp2})
            elif tmp2 == "None":
                locus_dict.update({i: "Orphan " + tmp1})
        if any(re.search("No_contig", tmp) for tmp in [tmp1, tmp2]):
            locus_dict.update({i: "No_contig"})
    result = pd.Series(locus_dict)
    return result


def lib_type(lib: str):
    """Dictionary of acceptable loci for library type."""
    librarydict = {
        "tr-ab": ["TRA", "TRB"],
        "tr-gd": ["TRG", "TRD"],
        "ig": ["IGH", "IGK", "IGL"],
    }
    return librarydict[lib]


def movecol(
    df: pd.DataFrame,
    cols_to_move: list = [],
    ref_col: str = "",
) -> pd.DataFrame:
    """A way to order columns."""
    # https://towardsdatascience.com/reordering-pandas-dataframe-columns-thumbs-down-on-standard-solutions-1ff0bc2941d5
    cols = df.columns.tolist()
    seg1 = cols[: list(cols).index(ref_col) + 1]
    seg2 = cols_to_move

    seg1 = [i for i in seg1 if i not in seg2]
    seg3 = [i for i in cols if i not in seg1 + seg2]
    return df[seg1 + seg2 + seg3]


def format_chain_status(locus_status):
    """Format chain status from locus status."""
    chain_status = []
    for ls in locus_status:
        if ("Orphan" in ls) and (re.search("TRB|IGH|TRD|VDJ", ls)):
            if not re.search("exception", ls):
                chain_status.append("Orphan VDJ")
            else:
                chain_status.append("Orphan VDJ-exception")
        elif ("Orphan" in ls) and (re.search("TRA|TRG|IGK|IGL|VJ", ls)):
            if not re.search("exception", ls):
                chain_status.append("Orphan VJ")
            else:
                chain_status.append("Orphan VJ-exception")
        elif re.search("exception|IgM/IgD", ls):
            chain_status.append("Extra pair-exception")
        elif re.search("Extra", ls):
            chain_status.append("Extra pair")
        elif re.search("ambiguous|None", ls):
            chain_status.append("ambiguous")
        else:
            chain_status.append("Single pair")
    return chain_status


def set_germline_env(
    germline: str | None = None,
    org: Literal["human", "mouse"] = "human",
    input_file: Path | str | None = None,
    db: Literal["imgt", "ogrdb"] = "imgt",
) -> tuple[dict[str, str], Path, Path]:
    """
    Set the paths to germline database and environment variables and relevant input files.

    Parameters
    ----------
    germline : str | None, optional
        path to germline database. None defaults to environmental variable $GERMLINE.
    org : Literal["human", "mouse"], optional
        organism for germline sequences.
    input_file : Path | str | None, optional
        path to input file.
    db : Literal["imgt", "ogrdb"], optional
        database to use. Defaults to imgt.
    Returns
    -------
    tuple[dict[str, str], Path, Path]
        environment dictionary and path to germline database.

    Raises
    ------
    KeyError
        if $GERMLINE environmental variable is not set.
    """
    env = os.environ.copy()
    if germline is None:
        try:
            gml = Path(env["GERMLINE"])
        except KeyError:
            raise KeyError(
                "Environmental variable $GERMLINE is missing. "
                "Please 'export GERMLINE=/path/to/database/germlines/'"
            )
        gml = gml / db / org / "vdj"
    else:
        gml = env["GERMLINE"] = Path(germline)
    if input_file is not None:
        input_file = Path(input_file)
    return env, gml, input_file


def set_igblast_env(
    igblast_db: Path | str | None = None,
    input_file: Path | str | None = None,
) -> tuple[dict[str, str], Path, Path]:
    """
    Set the igblast database and environment variables and relevant input files.

    Parameters
    ----------
    igblast_db : str | None, optional
        path to igblast database. None defaults to environmental variable $IGDATA.
    input_file : Path | str | None, optional
        path to input file.

    Returns
    -------
    tuple[dict[str, str], Path, Path]
        environment dictionary and path to igblast database.

    Raises
    ------
    KeyError
        if $IGDATA environmental variable is not set.
    """
    env = os.environ.copy()
    if igblast_db is None:
        try:
            igdb = Path(env["IGDATA"])
        except KeyError:
            raise KeyError(
                "Environmental variable $IGDATA is missing. "
                "Please 'export IGDATA=/path/to/database/igblast/'"
            )
    else:
        igdb = env["IGDATA"] = Path(igblast_db)
    if input_file is not None:
        input_file = Path(input_file)
    return env, igdb, input_file


def set_blast_env(
    blast_db: str | None = None,
    input_file: Path | str | None = None,
) -> tuple[dict[str, str], Path, Path]:
    """
    Set the blast database and environment variables and relevant input files.

    Parameters
    ----------
    blast_db : str | None, optional
        path to blast database. None defaults to environmental variable $BLASTDB.
    input_file : Path | str | None, optional
        path to input file.
    Returns
    -------
    tuple[dict[str, str], Path, Path]
        environment dictionary and path to igblast database.

    Raises
    ------
    KeyError
        if $BLASTDB environmental variable is not set.
    """
    env = os.environ.copy()
    if blast_db is None:
        try:
            bdb = Path(env["BLASTDB"])
        except KeyError:
            raise KeyError(
                "Environmental variable $BLASTDB is missing. "
                "Please 'export BLASTDB=/path/to/database/blast/'"
            )
    else:
        bdb = env["BLASTDB"] = Path(blast_db)
    if input_file is not None:
        input_file = Path(input_file)
    return env, bdb, input_file


def sum_col(vals: list) -> float | int:
    """Sum columns if not NaN."""
    if all(pd.isnull(vals)):
        return np.nan
    else:
        return sum(vals)


def check_data(
    data: list[Path | str] | Path | str, filename_prefix: list[str] | str | None
) -> tuple[list[str], list[str]]:
    """Quick check for data and filename prefixes"""
    if type(data) is not list:
        data = [data]
    if not isinstance(filename_prefix, list):
        filename_prefix = [filename_prefix]
        if len(filename_prefix) == 1:
            if len(data) > 1:
                filename_prefix = filename_prefix * len(data)
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for d in data]
    return data, filename_prefix


def check_same_celltype(clone_def1: str, clone_def2: str) -> bool:
    """Check if the first key is the same."""
    return clone_def1.split("_", 1)[0] == clone_def2.split("_", 1)[0]


def clear_h5file(filename: Path | str) -> None:
    """Little hack to overwrite an existing h5 file."""
    with h5py.File(filename, "w") as hf:
        for datasetname in hf.keys():
            del hf[datasetname]


def write_fasta(
    fasta_dict: dict[str, str], out_fasta: Path | str, overwrite=True
) -> None:
    """
    Generic fasta writer using fasta_iterator

    Parameters
    ----------
    fasta_dict : dict[str, str]
        dictionary containing fasta headers and sequences as keys and records respectively.
    out_fasta : str
        path to write fasta file to.
    overwrite : bool, optional
        whether or not to overwrite the output file (out_fasta).
    """
    if overwrite:
        fh = open(out_fasta, "w")
        fh.close()
    out = ""
    for l in fasta_dict:
        out = ">" + l + "\n" + fasta_dict[l] + "\n"
        write_output(out, out_fasta)


def write_output(out: str, file: Path | str) -> None:
    """General line writer."""
    fh = open(file, "a")
    fh.write(out)
    fh.close()
