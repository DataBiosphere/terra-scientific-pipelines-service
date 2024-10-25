# tests/commands/pipelines_tests.py

from unittest.mock import patch
from click.testing import CliRunner
from teaspoons.commands.pipelines_commands import pipelines
from teaspoons_client.models.pipeline import Pipeline
from teaspoons_client.models.pipeline_with_details import PipelineWithDetails
from teaspoons_client.exceptions import ApiException


def test_list_pipelines():
    runner = CliRunner()

    with patch(
        "teaspoons.commands.pipelines_commands.pipelines_logic.list_pipelines",
        return_value=[Pipeline(pipeline_name="test_pipeline_1", display_name="test_display_name_1", description="test_description_1"), 
                      Pipeline(pipeline_name="test_pipeline_2", display_name="test_display_name_2", description="test_description_2")]
    ) as mock_list_pipelines:
        result = runner.invoke(pipelines, ["list"])

        # Assert the command executed successfully
        assert result.exit_code == 0
        # Assert the logic function was called
        mock_list_pipelines.assert_called_once()
        # Assert the output is formatted correctly
        assert "Found 2 available pipelines:\n\ttest_pipeline_1 - test_description_1\n\ttest_pipeline_2 - test_description_2" in result.output


def test_get_info_success():
    runner = CliRunner()

    with patch(
        "teaspoons.commands.pipelines_commands.pipelines_logic.get_pipeline_info",
        return_value=PipelineWithDetails(pipeline_name="test_pipeline", description="test_description", display_name="test_display_name", type="test_type", inputs=[])
    ) as mock_get_pipeline_info:
        result = runner.invoke(pipelines, ["get-info", "test_pipeline"])

        # Assert the command executed successfully
        assert result.exit_code == 0
        # Assert the logic function was called with the correct argument
        mock_get_pipeline_info.assert_called_once_with("test_pipeline")
        # Assert the output contains the pipeline name (could test other things here)
        assert "test_pipeline" in result.output


def test_get_info_missing_argument():
    runner = CliRunner()

    # Assert the command raises a PipelineApi exception
    result = runner.invoke(pipelines, ["get-info"])

    # Assert the command failed due to missing argument
    assert result.exit_code != 0
    assert "Error: Missing argument 'PIPELINE_NAME'" in result.output


def test_get_info_bad_pipeline_name():
    runner = CliRunner()

    with patch("teaspoons.commands.pipelines_commands.pipelines_logic.get_pipeline_info", side_effect=ApiException(404, "Pipeline not found")):
        result = runner.invoke(pipelines, ["get-info", "bad_pipeline_name"])

        # Assert the command failed due to an ApiException and the pipeline not being found
        assert result.exit_code != 0
        assert isinstance(result.exception, ApiException)
        assert "Pipeline not found" in result.exception.reason
