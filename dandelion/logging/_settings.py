import inspect
import sys
from contextlib import contextmanager
from enum import IntEnum
from pathlib import Path
from time import time
from logging import getLevelName
from typing import Any, Union, Optional, Iterable, TextIO
from typing import Tuple, List, ContextManager

from . import logger as logging
from .logger import _set_log_level, _set_log_file, _RootLogger
from ..utilities._utilities import Literal

_VERBOSITY_TO_LOGLEVEL = {
    'error': 'ERROR',
    'warning': 'WARNING',
    'info': 'INFO',
    'hint': 'HINT',
    'debug': 'DEBUG',
}
# Python 3.7 ensures iteration order
for v, level in enumerate(list(_VERBOSITY_TO_LOGLEVEL.values())):
    _VERBOSITY_TO_LOGLEVEL[v] = level


class Verbosity(IntEnum):
    error = 0
    warn = 1
    info = 2
    hint = 3
    debug = 4

    @property
    def level(self) -> int:
        # getLevelName(str) returns the int level…
        return getLevelName(_VERBOSITY_TO_LOGLEVEL[self])

    @contextmanager
    def override(self, verbosity: "Verbosity") -> ContextManager["Verbosity"]:
        """\
        Temporarily override verbosity
        """
        settings.verbosity = verbosity
        yield self
        settings.verbosity = self


def _type_check(var: Any, varname: str, types: Union[type, Tuple[type, ...]]):
    if isinstance(var, types):
        return
    if isinstance(types, type):
        possible_types_str = types.__name__
    else:
        type_names = [t.__name__ for t in types]
        possible_types_str = "{} or {}".format(", ".join(type_names[:-1]),
                                               type_names[-1])
    raise TypeError(f"{varname} must be of type {possible_types_str}")


class Config:
    """
    Config manager
    """

    def __init__(
        self,
        *,
        verbosity: str = "warning",
        logfile: Union[str, Path, None] = None,
    ):
        # logging
        self._root_logger = _RootLogger(logging.INFO)  # level will be replaced
        self.logfile = logfile
        self.verbosity = verbosity
        self._start = time()
        """Time when the settings module is first imported."""

        self._previous_time = self._start
        """Variable for timing program parts."""

    @property
    def verbosity(self) -> Verbosity:
        """
        Verbosity level (default `warning`)
        Level 0: only show 'error' messages.
        Level 1: also show 'warning' messages.
        Level 2: also show 'info' messages.
        Level 3: also show 'hint' messages.
        Level 4: also show very detailed progress for 'debug'ging.
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: Union[Verbosity, int, str]):
        verbosity_str_options = [
            v for v in _VERBOSITY_TO_LOGLEVEL if isinstance(v, str)
        ]
        if isinstance(verbosity, Verbosity):
            self._verbosity = verbosity
        elif isinstance(verbosity, int):
            self._verbosity = Verbosity(verbosity)
        elif isinstance(verbosity, str):
            verbosity = verbosity.lower()
            if verbosity not in verbosity_str_options:
                raise ValueError(
                    f"Cannot set verbosity to {verbosity}. "
                    f"Accepted string values are: {verbosity_str_options}")
            else:
                self._verbosity = Verbosity(
                    verbosity_str_options.index(verbosity))
        else:
            _type_check(verbosity, "verbosity", (str, int))
        _set_log_level(self, _VERBOSITY_TO_LOGLEVEL[self._verbosity])

    @property
    def logpath(self) -> Optional[Path]:
        """\
        The file path `logfile` was set to.
        """
        return self._logpath

    @logpath.setter
    def logpath(self, logpath: Union[str, Path, None]):
        _type_check(logpath, "logfile", (str, Path))
        # set via “file object” branch of logfile.setter
        self.logfile = Path(logpath).open('a')
        self._logpath = Path(logpath)

    @property
    def logfile(self) -> TextIO:
        """\
        The open file to write logs to.
        Set it to a :class:`~pathlib.Path` or :class:`str` to open a new one.
        The default `None` corresponds to :obj:`sys.stdout` in jupyter notebooks
        and to :obj:`sys.stderr` otherwise.
        For backwards compatibility, setting it to `''` behaves like setting it to `None`.
        """
        return self._logfile

    @logfile.setter
    def logfile(self, logfile: Union[str, Path, TextIO, None]):
        if not hasattr(logfile, 'write') and logfile:
            self.logpath = logfile
        else:  # file object
            if not logfile:  # None or ''
                logfile = sys.stdout if self._is_run_from_ipython(
                ) else sys.stderr
            self._logfile = logfile
            self._logpath = None
            _set_log_file(self)

    @staticmethod
    def _is_run_from_ipython():
        """Determines whether we're currently in IPython."""
        import builtins

        return getattr(builtins, "__IPYTHON__", False)

    def __str__(self) -> str:
        return '\n'.join(f'{k} = {v!r}' for k, v in inspect.getmembers(self)
                         if not k.startswith("_") and not k == 'getdoc')


settings = Config()
