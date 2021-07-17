#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2021-04-03 16:46:13
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-07-17 21:32:49

import traceback
from pathlib import Path

here = Path(__file__).parent.parent

try:
    from setuptools_scm import get_version
    import pytoml

    proj = pytoml.loads((here.parent / "pyproject.toml").read_text())
    metadata = proj["tool"]["flit"]["metadata"]

    # __version__ = get_version(root="../..", relative_to=__file__, **proj["tool"]["setuptools_scm"]).split("+")[0]
    __version__ = metadata["version"]
    __author__ = metadata["author"]
    __email__ = metadata["author-email"]
    __url__ = metadata["home-page"]
    __docs = metadata["documentation"]

except (ImportError, LookupError, FileNotFoundError):
    try:
        from importlib.metadata import metadata as m
    except ImportError:  # < Python 3.8: Use backport module
        from importlib_metadata import metadata as m

    metadata = m(here.name)
    __version__ = metadata["Version"]
    __author__ = metadata["Author"]
    __email__ = metadata["Author-email"]
    __url__ = metadata["Home-page"]
    __docs__ = metadata["Documentation"]


def within_flit():
    for frame in traceback.extract_stack():
        if frame.name == "get_docstring_and_version_via_import":
            return True
    return False
