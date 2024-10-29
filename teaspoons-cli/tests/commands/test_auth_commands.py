# tests/commands/test_auth_commands.py

from click.testing import CliRunner
from mockito import when, verify
from teaspoons.commands import auth_commands


def test_login():
    runner = CliRunner()

    when(auth_commands.auth_logic).check_local_token_and_fetch_if_needed()

    result = runner.invoke(auth_commands.auth, ["login"])

    assert result.exit_code == 0
    verify(auth_commands.auth_logic).check_local_token_and_fetch_if_needed()


def test_logout():
    runner = CliRunner()

    when(auth_commands.auth_logic).clear_local_token()

    result = runner.invoke(auth_commands.auth, ["logout"])

    assert result.exit_code == 0
    verify(auth_commands.auth_logic).clear_local_token()
