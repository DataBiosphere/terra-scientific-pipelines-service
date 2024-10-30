# tests/test_cli.py

from click.testing import CliRunner

from teaspoons import cli


def test_cli():
    runner = CliRunner()

    result = runner.invoke(cli.cli, ["-h"])

    # expect the following packages
    expected_commands = [
        "auth",
        "pipelines"
    ]
    for command in expected_commands:
        assert command in result.output
