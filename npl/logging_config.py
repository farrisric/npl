"""User-facing logging setup for npl.

Library modules only call ``logging.getLogger(__name__)``; importing npl
attaches no handlers. Applications that want log output call
:func:`setup_logging` once at startup.
"""
import logging


def setup_logging(level=logging.INFO, file=None, stream=True):
    """Attach handlers to the top-level ``npl`` logger, never the root logger.

    Args:
        level: Log level for the npl logger.
        file: Optional path to a log file, e.g. ``'npl.log'``.
        stream: If True, also emit to stderr.

    Returns:
        The configured ``npl`` logger.
    """
    logger = logging.getLogger('npl')
    logger.setLevel(level)
    logger.propagate = False

    formatter = logging.Formatter('%(levelname)s - %(message)s')

    # Drop handlers a previous call attached, so calling again reconfigures.
    for h in list(logger.handlers):
        if getattr(h, '_npl_managed', False):
            logger.removeHandler(h)
            h.close()

    handlers = []
    if stream:
        handlers.append(logging.StreamHandler())
    if file is not None:
        handlers.append(logging.FileHandler(file))
    for h in handlers:
        h.setFormatter(formatter)
        h._npl_managed = True
        logger.addHandler(h)

    return logger
