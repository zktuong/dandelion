#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2021-04-03 16:46:13
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-07-18 17:29:33

# try:
#   from .version import __version__
# except ImportError:
#   __version__ = ''
from pathlib import Path

here = Path(__file__).parent

try:
    from setuptools_scm import get_version
    import pytoml
    proj = pytoml.loads((here.parent / "pyproject.toml").read_text())

    __version__ = get_version(root="..", relative_to=__file__, **proj["tool"]["setuptools_scm"])
    __author__ = "Zewen Kelvin Tuong"
    __email__ = "kt16@sanger.ac.uk"
    __classifiers__ = [
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: R",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
    ]
except:
    try:
        from importlib.metadata import metadata as m
    except ImportError:
        from importlib_metadata import metadata as m

    metadata = m(here.name)
    __version__ = metadata["Version"]
    __author__ = metadata["Author"]
    __email__ = metadata["Author-email"]
