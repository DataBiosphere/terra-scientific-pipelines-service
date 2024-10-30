# log.py

import logging


def configure_logging(debug: bool):
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=log_level, format="%(message)s")
