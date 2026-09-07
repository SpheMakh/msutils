"""Console logging for the two console scripts.

The library half of msutils never touches a handler: a module takes its
logger with `logging.getLogger(__name__)` and emits. Attaching the handler is
the application's job, so `configure_logging` is called from `msutils.cli`
and `msutils.gains.cli` and from nowhere else. `import msutils` therefore
configures nothing, and a caller embedding msutils keeps whatever logging
they have already set up.

This used to be `create_logger`, called at import by fourteen library modules,
so importing any one of them attached a StreamHandler to the process.
"""

import logging

from ._package import LOGGER


def configure_logging(level=logging.INFO):
    """Attach msutils' console handler, at `level`.

    Idempotent: a second call re-levels the handler already attached rather
    than stacking another one, which matters when both console scripts run in
    one process (a test, or a cab importing both).

    Args:
        level: A logging level, or a level name.

    Returns:
        logging.Logger: The `msutils` logger, configured.
    """
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    log = logging.getLogger(LOGGER)
    log.setLevel(level)

    if not any(not isinstance(handler, logging.NullHandler) for handler in log.handlers):
        console = logging.StreamHandler()
        console.setFormatter(
            logging.Formatter("%(name)s - %(asctime)s %(levelname)s - %(message)s")
        )
        log.addHandler(console)
    for handler in log.handlers:
        handler.setLevel(level)

    return log
