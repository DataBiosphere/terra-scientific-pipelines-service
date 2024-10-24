# tests/commands/test_auth_commands.py

from unittest.mock import patch
from click.testing import CliRunner
from teaspoons.commands.auth_commands import auth


def test_login():
    runner = CliRunner()

    with patch(
        "teaspoons.commands.auth_commands.auth_logic.check_local_token_and_fetch_if_needed"
    ) as mock_check_local_token_and_fetch_if_needed:
        result = runner.invoke(auth, ["login"])

        assert result.exit_code == 0
        mock_check_local_token_and_fetch_if_needed.assert_called_once()


def test_logout():
    runner = CliRunner()
    with patch(
        "teaspoons.commands.auth_commands.auth_logic.clear_local_token"
    ) as mock_clear_local_token:
        result = runner.invoke(auth, ["logout"])

        assert result.exit_code == 0
        mock_clear_local_token.assert_called_once()
