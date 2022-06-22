#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2022-06-21 21:40:09
"""utilities module."""
import numpy as np
import os
import pandas as pd
import re
import warnings

from airr import RearrangementSchema
from collections import defaultdict
from subprocess import run
from typing import Sequence, Tuple, Dict, Union, Optional

# for compatibility with python>=3.10
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable


try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            """LiteralMeta class."""

            def __getitem__(self, values):
                """Return Literal."""
                if not isinstance(values, tuple):
                    values = (values,)
                return type("Literal_", (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            """Literal type."""

            pass


class Tree(defaultdict):
    """Create a recursive defaultdict."""

    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value


def dict_from_table(meta: pd.DataFrame, columns: Tuple[str, str]) -> Dict:
    """
    Generate a dictionary from a dataframe.

    Parameters
    ----------
    meta
        pandas dataframe or file path
    columns
        column names in dataframe

    Returns
    -------
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


def clean_nan_dict(d: Dict) -> Dict:
    """
    Remove nan from dictionary.

    Parameters
    ----------
    d : Dict
        dictionary

    Returns
    -------
        dictionary with no NAs.
    """
    return {k: v for k, v in d.items() if v is not np.nan}


def flatten(l: Sequence) -> Sequence:
    """
    Flatten a list-in-list-in-list.

    Parameters
    ----------
    l : Sequence
        a list-in-list list

    Returns
    -------
        a flattened list.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def makeblastdb(ref: str):
    """
    Run makeblastdb on constant region fasta file.

    Wrapper for makeblastdb.

    Parameters
    ----------
    ref : str
        constant region fasta file
    """
    cmd = ["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", ref]
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


def is_categorical(array_like) -> bool:
    """Check if categorical."""
    return array_like.dtype.name == "category"


def type_check(dataframe, key) -> bool:
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
    s,
    filename_prefix: Optional[str] = None,
    endswith: Optional[str] = None,
    subdir: Optional[str] = None,
):
    """Check filepath."""
    if filename_prefix is None:
        filename_pre = "filtered"
    else:
        filename_pre = filename_prefix

    if endswith is None:
        ends_with = ""
    else:
        ends_with = endswith

    filePath = None
    if os.path.isfile(str(s)) and str(s).endswith(ends_with):
        filePath = s
    elif os.path.isdir(str(s)):
        files = os.listdir(s)
        for file in files:
            out_ = s.rstrip("/") + "/"
            if os.path.isdir(out_ + os.path.basename(file)):
                if file == "dandelion":
                    if subdir is None:
                        out_ = out_ + os.path.basename(file) + "/"
                    else:
                        out_ = (
                            out_ + os.path.basename(file) + "/" + subdir + "/"
                        )
                    for x in os.listdir(out_):
                        if x.endswith(ends_with):
                            if (
                                str(x).split(ends_with)[0]
                                == filename_pre + "_contig"
                            ):
                                filePath = out_ + x
                            else:
                                continue
            else:
                continue
    return filePath


def check_fastapath(fasta, filename_prefix: Optional[str] = None):
    """Check fastapath."""
    if filename_prefix is None:
        filename_pre = "filtered"
    else:
        filename_pre = filename_prefix

    filePath = None
    if os.path.isfile(str(fasta)) and str(fasta).endswith(".fasta"):
        filePath = fasta
    elif os.path.isdir(str(fasta)):
        files = os.listdir(fasta)
        for file in files:
            out_ = fasta.rstrip("/") + "/"
            if str(file).endswith(".fasta"):
                if str(file).split(".fasta")[0] == filename_pre + "_contig":
                    filePath = out_ + os.path.basename(file)
                else:
                    continue
    return filePath


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

        def __init__(self, obj, *args):
            self.obj = obj

        def __lt__(self, other):
            """Less than."""
            return mycmp(self.obj, other.obj) < 0

        def __gt__(self, other):
            """Greater than."""
            return mycmp(self.obj, other.obj) > 0

        def __eq__(self, other):
            """Equal."""
            return mycmp(self.obj, other.obj) == 0

        def __le__(self, other):
            """Less than or equal."""
            return mycmp(self.obj, other.obj) <= 0

        def __ge__(self, other):
            """Greater than or equal."""
            return mycmp(self.obj, other.obj) >= 0

        def __ne__(self, other):
            """Not equal."""
            return mycmp(self.obj, other.obj) != 0

    return K


def not_same_call(a, b, pattern):
    """Utility function to check if a == b in terms of pattern."""
    return (re.search(pattern, a) and not re.search(pattern, b)) or (
        re.search(pattern, b) and not re.search(pattern, a)
    )


def same_call(a, b, c, pattern):
    """Utility function to check if a == b == c in terms of pattern."""
    queries = [a, b, c]
    queries = [q for q in queries if pd.notnull(q)]
    return all([re.search(pattern, x) for x in queries])


def present(x):
    """Utility function to check if x is not null or blank."""
    return pd.notnull(x) and x != ""


def check_missing(x):
    """Utility function to check if x is null or blank."""
    return pd.isnull(x) or x == ""


def all_missing(x):
    """Utility function to check if all x is not null or blank."""
    return all(pd.isnull(x)) or all(x == "")


def return_mix_dtype(data):
    """Utility function to return mixed dtypes columns."""
    check = [
        c
        for c in data.columns
        if pd.api.types.infer_dtype(data[c]).startswith("mixed")
    ]
    return check


def sanitize_data(data, ignore="clone_id"):
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
                data[d].replace(
                    [None, np.nan, pd.NA, "nan", ""], "", inplace=True
                )
                if RearrangementSchema.properties[d]["type"] == "integer":
                    data[d] = [
                        int(x) if present(x) else ""
                        for x in pd.to_numeric(data[d])
                    ]
            else:
                data[d].replace(
                    [None, pd.NA, np.nan, "nan", ""], np.nan, inplace=True
                )
        else:
            if d != ignore:
                try:
                    data[d] = pd.to_numeric(data[d])
                except:
                    data[d].replace(
                        to_replace=[None, np.nan, pd.NA, "nan", ""],
                        value="",
                        inplace=True,
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
        pd.Series(["duplicate_count", "productive"]).isin(data.columns).all()
    ):  # sort so that the productive contig with the largest umi is first
        data.sort_values(
            by=["productive", "duplicate_count"], inplace=True, ascending=False
        )

    # check if airr-standards is happy
    validate_airr(data)
    return data


def sanitize_blastn(data):
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
                data[d].replace(
                    [None, np.nan, pd.NA, "nan", ""], "", inplace=True
                )
                if RearrangementSchema.properties[d]["type"] == "integer":
                    data[d] = [
                        int(x) if present(x) else ""
                        for x in pd.to_numeric(data[d])
                    ]
            else:
                data[d].replace(
                    [None, pd.NA, np.nan, "nan", ""], np.nan, inplace=True
                )
        else:
            try:
                data[d] = pd.to_numeric(data[d])
            except:
                data[d].replace(
                    to_replace=[None, np.nan, pd.NA, "nan", ""],
                    value="",
                    inplace=True,
                )
    return data


def sanitize_data_for_saving(data):
    """Quick sanitize dtypes for saving."""
    tmp = data.copy()
    for d in tmp:
        if d in RearrangementSchema.properties:
            if RearrangementSchema.properties[d]["type"] in [
                "string",
                "boolean",
            ]:
                tmp[d].replace(
                    [None, np.nan, pd.NA, "nan", ""], "", inplace=True
                )
            if RearrangementSchema.properties[d]["type"] in [
                "integer",
                "number",
            ]:
                tmp[d].replace(
                    [None, np.nan, pd.NA, "nan", ""], np.nan, inplace=True
                )
        else:
            try:
                tmp[d] = pd.to_numeric(tmp[d])
            except:
                tmp[d].replace(
                    [None, pd.NA, np.nan, "nan", ""], "", inplace=True
                )
    return tmp


def validate_airr(data):
    """Validate dtypes in airr table."""
    tmp = data.copy()
    int_columns = []
    for d in tmp:
        try:
            tmp[d].replace(np.nan, pd.NA).astype("Int64")
            int_columns.append(d)
        except:
            pass
    bool_columns = [
        "rev_comp",
        "productive",
        "vj_in_frame",
        "stop_codon",
        "complete_vdj",
    ]
    str_columns = list(tmp.dtypes[tmp.dtypes == "object"].index)
    columns = [
        c
        for c in list(set(int_columns + str_columns + bool_columns))
        if c in tmp
    ]
    if len(columns) > 0:
        for c in columns:
            tmp[c].fillna("", inplace=True)
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


def check_travdv(data):
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


def load_data(obj: Optional[Union[pd.DataFrame, str]]) -> pd.DataFrame:
    """
    Read in or copy dataframe object and set sequence_id as index without dropping.

    Parameters
    ----------
    obj : DataFrame, str
        file path to .tsv file or pandas DataFrame object.

    Returns
    -------
    pandas DataFrame object.
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
            obj_.set_index("sequence_id", drop=False, inplace=True)
            if "cell_id" not in obj_.columns:
                obj_["cell_id"] = [
                    c.split("_contig")[0] for c in obj_["sequence_id"]
                ]
        else:
            raise KeyError("'sequence_id' not found in columns of input")

        if "umi_count" in obj_.columns:
            if "duplicate_count" not in obj_.columns:
                obj_.rename(
                    columns={"umi_count": "duplicate_count"}, inplace=True
                )

        return obj_


class ContigDict(dict):
    """Class Object to extract the contigs as a dictionary."""

    def __setitem__(self, key, value):
        """Standard __setitem__."""
        super().__setitem__(key, value)

    def __hash__(self):
        """Make it hashable."""
        return hash(tuple(self))


class Contig:
    """Class Object to hold contig."""

    def __init__(self, contig, mapper=None):
        if mapper is not None:
            mapper.update({k: k for k in contig.keys() if k not in mapper})
            self._contig = ContigDict(
                {mapper[key]: vals for (key, vals) in contig.items()}
            )
        else:
            self._contig = ContigDict(contig)

    @property
    def contig(self):
        """Contig slot."""
        return self._contig


def mask_dj(data, filename_prefix, d_evalue_threshold, j_evalue_threshold):
    """Mask d/j assignment."""
    for i in range(0, len(data)):
        filePath = check_filepath(
            data[i],
            filename_prefix=filename_prefix[i],
            endswith="_igblast_db-pass.tsv",
        )
        if filePath is None:
            raise FileNotFoundError(
                "Path to .tsv file for {} is unknown. ".format(data[i])
                + "Please specify path to reannotated .tsv file or folder containing reannotated .tsv file."
            )

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


def write_airr(data, save):
    """Save as airr formatted file."""
    data = sanitize_data(data)
    data.to_csv(save, sep="\t", index=False)


def write_blastn(data, save):
    """Write blast output."""
    data = sanitize_blastn(data)
    data.to_csv(save, sep="\t", index=False)


## from skbio==0.5.6
def _validate_counts_vector(counts, suppress_cast=False):
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


def deprecated(details, deprecated_in, removed_in):
    """Decorator to mark a function as deprecated"""

    def deprecated_decorator(func):
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
    suffix_h: str = "_VDJ",
    suffix_l: str = "_VJ",
) -> list:
    """Extract v/d/j/c call values from data."""
    call_dict = {
        "Multi": "Multi",
        "None": "None",
        "": "None",
        "unassigned": "None",
    }
    if suffix_l is not None:
        call_1 = {
            x[0]: x[1] if present(x[1]) else "None"
            for x, y in zip(
                metadata[call + suffix_h].items(),
                list(metadata[call + suffix_l]),
            )
        }
        call_2 = {
            x[0]: x[1] if present(x[1]) else "None"
            for x, y in zip(
                metadata[call + suffix_l].items(),
                list(metadata[call + suffix_h]),
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
            for x, y in metadata[call + suffix_h].items()
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
    metadata: pd.DataFrame, suffix_h: str = "_VDJ", suffix_l: str = "_VJ"
) -> list:
    """Extract locus call value from data."""
    locus_1 = {
        x[0]: x[1] if present(x[1]) else y
        for x, y in zip(
            metadata["locus" + suffix_h].items(),
            list(metadata["locus" + suffix_l]),
        )
    }
    locus_2 = {
        x[0]: x[1] if present(x[1]) else y
        for x, y in zip(
            metadata["locus" + suffix_l].items(),
            list(metadata["locus" + suffix_h]),
        )
    }
    multi_1 = {
        x: "Multi" for x, y in metadata["locus" + suffix_h].items() if "|" in y
    }
    multi_2 = {
        x: "Multi" for x, y in metadata["locus" + suffix_l].items() if "|" in y
    }
    locus_1.update(multi_1)
    locus_2.update(multi_2)
    result = [
        str(x) + " + " + str(y) if str(x) != str(y) else str(x) + "_only"
        for x, y in zip(locus_1.values(), locus_2.values())
    ]
    result = [x if "Multi" not in x else "Multi" for x in result]
    return result


def format_productive(
    metadata: pd.DataFrame, suffix_h: str = "_VDJ", suffix_l: str = "_VJ"
) -> list:
    """Extract productive value from data."""
    productive_1 = {
        x[0]: x[1] if present(x[1]) else "None"
        for x, y in zip(
            metadata["productive" + suffix_h].items(),
            list(metadata["productive" + suffix_l]),
        )
    }
    productive_2 = {
        x[0]: x[1] if present(x[1]) else "None"
        for x, y in zip(
            metadata["productive" + suffix_l].items(),
            list(metadata["productive" + suffix_h]),
        )
    }
    multi_1 = {
        x: "Multi"
        for x, y in metadata["productive" + suffix_h].items()
        if "|" in y
    }
    multi_2 = {
        x: "Multi"
        for x, y in metadata["productive" + suffix_l].items()
        if "|" in y
    }
    productive_1.update(multi_1)
    productive_2.update(multi_2)
    result = [
        str(x) + " + " + str(y)
        for x, y in zip(productive_1.values(), productive_2.values())
    ]
    # result = [x if 'Multi' not in x else 'Multi' for x in result]
    return result


def sum_col(vals):
    """Sum columns if not NaN."""
    if all(pd.isnull(vals)):
        return np.nan
    else:
        return sum(vals)
