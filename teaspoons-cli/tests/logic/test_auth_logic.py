# tests/logic/test_auth_logic.py

import pytest
from mockito import when, mock, verify

from teaspoons.logic import auth_logic


@pytest.fixture
def mock_cli_config(unstub):
    config = mock({"token_file": "mock_token_file", "client_info": "mock_client_info"})
    when(auth_logic).CliConfig(...).thenReturn(config)
    yield config
    unstub()


def test_check_local_token_and_fetch_if_needed_already_authenticated(
    mock_cli_config, unstub
):
    when(auth_logic)._load_local_token(...).thenReturn("valid_token")
    when(auth_logic)._validate_token(...).thenReturn(True)

    auth_logic.check_local_token_and_fetch_if_needed()

    # Verify the token was loaded and validated
    verify(auth_logic)._load_local_token("mock_token_file")
    verify(auth_logic)._validate_token("valid_token")


def test_check_local_token_and_fetch_if_needed_fetch_new_token(mock_cli_config, unstub):
    when(auth_logic)._load_local_token(...).thenReturn(None)
    when(auth_logic).get_access_token_with_browser_open(...).thenReturn("new_token")
    when(auth_logic)._save_local_token(...).thenReturn(None)

    auth_logic.check_local_token_and_fetch_if_needed()

    # Verify the sequence of operations
    verify(auth_logic)._load_local_token("mock_token_file")
    verify(auth_logic).get_access_token_with_browser_open("mock_client_info")
    verify(auth_logic)._save_local_token("mock_token_file", "new_token")


def test_clear_local_token(mock_cli_config, unstub):
    when(auth_logic)._clear_local_token(...).thenReturn(None)

    auth_logic.clear_local_token()

    # Verify the token was cleared and message was printed
    verify(auth_logic)._clear_local_token("mock_token_file")
