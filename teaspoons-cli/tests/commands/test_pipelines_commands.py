# tests/commands/pipelines_tests.py

from unittest.mock import patch
from click.testing import CliRunner
from teaspoons.commands.pipelines_commands import pipelines


def test_list_pipelines():
    runner = CliRunner()

    # Mock the logic function as it's called by the command to bypass the decorator as well
    with patch(
        "teaspoons.commands.pipelines_commands.pipelines_logic.list_pipelines"
    ) as mock_list_pipelines:
        result = runner.invoke(pipelines, ["list"])

        print(result.output)
        # Assert the command executed successfully
        assert result.exit_code == 0
        # Assert the logic function was called
        mock_list_pipelines.assert_called_once()


def test_get_info_success():
    runner = CliRunner()

    with patch(
        "teaspoons.commands.pipelines_commands.pipelines_logic.get_pipeline_info"
    ) as mock_get_pipeline_info:
        result = runner.invoke(pipelines, ["get-info", "test_pipeline"])

        # Assert the command executed successfully
        assert result.exit_code == 0
        # Assert the logic function was called with the correct argument
        mock_get_pipeline_info.assert_called_once_with("test_pipeline")


def test_get_info_missing_argument():
    runner = CliRunner()

    result = runner.invoke(pipelines, ["get-info"])

    # Assert the command failed due to missing argument
    assert result.exit_code != 0
    assert "Error: Missing argument 'PIPELINE_NAME'" in result.output
