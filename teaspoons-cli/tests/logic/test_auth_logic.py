# tests/logic/test_auth_logic.py

import pytest

from unittest.mock import patch
from teaspoons.logic.auth_logic import (
    check_local_token_and_fetch_if_needed,
    clear_local_token,
)


@pytest.fixture
def mock_cli_config():
    with patch("teaspoons.logic.auth_logic.CliConfig") as MockCliConfig:
        mock_config = MockCliConfig.return_value
        mock_config.token_file = "mock_token_file"
        mock_config.client_info = "mock_client_info"
        yield mock_config


def test_check_local_token_and_fetch_if_needed_already_authenticated(mock_cli_config):
    with patch(
        "teaspoons.logic.auth_logic._load_local_token", return_value="valid_token"
    ) as mock_load_token, patch(
        "teaspoons.logic.auth_logic._validate_token", return_value=True
    ) as mock_validate_token, patch(
        "builtins.print"
    ) as mock_print:

        check_local_token_and_fetch_if_needed()

        mock_load_token.assert_called_once_with("mock_token_file")
        mock_validate_token.assert_called_once_with("valid_token")
        mock_print.assert_called_once_with("Already authenticated")


def test_check_local_token_and_fetch_if_needed_fetch_new_token(mock_cli_config):
    with patch(
        "teaspoons.logic.auth_logic._load_local_token", return_value=None
    ), patch(
        "teaspoons.logic.auth_logic.get_access_token_with_browser_open",
        return_value="new_token",
    ), patch(
        "teaspoons.logic.auth_logic._save_local_token"
    ) as mock_save_token:

        check_local_token_and_fetch_if_needed()

        mock_save_token.assert_called_once_with("mock_token_file", "new_token")


def test_clear_local_token(mock_cli_config):
    with patch(
        "teaspoons.logic.auth_logic._clear_local_token"
    ) as mock_clear_token, patch("builtins.print") as mock_print:

        clear_local_token()

        mock_clear_token.assert_called_once_with("mock_token_file")
        mock_print.assert_called_once_with("Logged out")
