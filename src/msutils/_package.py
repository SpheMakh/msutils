"""Package-level constants, kept out of the package initializer.

An `__init__.py` that holds implementation is what RUF067
(`non-empty-init-module`) is about: importing the package then runs it. The
initializer re-exports what belongs on `msutils`; this is where the two
things that are not re-exports live.
"""

import logging
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _version

#: Read off the installed distribution rather than restated here, so it cannot
#: drift from `pyproject.toml`. The two console scripts take it from here too;
#: they each used to carry their own copy of this block.
try:
    __version__ = _version("msutils")
except PackageNotFoundError:  # running from a source tree, not installed
    __version__ = "0.0.0.dev0"

#: The root of the logger hierarchy every module emits through -- each takes
#: its own `logging.getLogger(__name__)`, which is a child of this.
LOGGER = "msutils"

# Library convention: modules only ever *emit*, and never attach a handler, so
# importing msutils configures nothing. The NullHandler is what makes that
# silence deliberate rather than leaving logging's last-resort stderr echo of
# unhandled WARNING+ records to decide. The console handler belongs to the
# application -- see `msutils._log.configure_logging`, called by both console
# scripts.
logging.getLogger(LOGGER).addHandler(logging.NullHandler())
