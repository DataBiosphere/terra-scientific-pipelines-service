# logic/auth_logic.py

import logging
from auth_helper import (
    get_access_token_with_browser_open,
    _validate_token,
    _save_local_token,
    _load_local_token,
    _clear_local_token,
)
from config import CliConfig

LOGGER = logging.getLogger(__name__)


def check_local_token_and_fetch_if_needed():
    """Authenticate with Teaspoons via browser login to Terra b2c"""
    cli_config = CliConfig()  # initialize the config from environment variables
    token = _load_local_token(cli_config.token_file)
    if token and _validate_token(token):
        LOGGER.info("Already authenticated")
        return
    token = get_access_token_with_browser_open(cli_config.client_info)
    _save_local_token(cli_config.token_file, token)


def clear_local_token():
    """Clear the local authentication token"""
    cli_config = CliConfig()  # initialize the config from environment variables
    _clear_local_token(cli_config.token_file)
    LOGGER.info("Logged out")
