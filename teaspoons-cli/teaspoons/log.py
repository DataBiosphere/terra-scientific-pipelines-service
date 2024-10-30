# log.py

import logging


def configure_logging(debug):
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=log_level, format="%(message)s")


