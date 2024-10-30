# tests/commands/pipelines_tests.py

import logging
from click.testing import CliRunner
from mockito import when, verify
from teaspoons.commands import pipelines_commands
from teaspoons_client.models.pipeline import Pipeline
from teaspoons_client.models.pipeline_with_details import PipelineWithDetails
from teaspoons_client.exceptions import ApiException


LOGGER = logging.getLogger(__name__)


def test_list_pipelines(caplog):
    runner = CliRunner()
    test_pipelines = [
        Pipeline(
            pipeline_name="test_pipeline_1",
            display_name="test_display_name_1",
            description="test_description_1",
        ),
        Pipeline(
            pipeline_name="test_pipeline_2",
            display_name="test_display_name_2",
            description="test_description_2",
        ),
    ]

    when(pipelines_commands.pipelines_logic).list_pipelines().thenReturn(test_pipelines)

    with caplog.at_level(logging.DEBUG):
        result = runner.invoke(pipelines_commands.pipelines, ["list"])

    assert result.exit_code == 0
    verify(pipelines_commands.pipelines_logic).list_pipelines()
    assert "Found 2 available pipelines:" in caplog.text
    assert "test_pipeline_1" in caplog.text
    assert "test_pipeline_2" in caplog.text


def test_get_info_success(caplog, unstub):
    runner = CliRunner()
    test_pipeline = PipelineWithDetails(
        pipeline_name="test_pipeline",
        description="test_description",
        display_name="test_display_name",
        type="test_type",
        inputs=[],
    )

    when(pipelines_commands.pipelines_logic).get_pipeline_info(
        "test_pipeline"
    ).thenReturn(test_pipeline)

    with caplog.at_level(logging.DEBUG):
        result = runner.invoke(
            pipelines_commands.pipelines, ["get-info", "test_pipeline"]
        )

    assert result.exit_code == 0
    verify(pipelines_commands.pipelines_logic).get_pipeline_info("test_pipeline")
    assert "test_pipeline" in caplog.text

    unstub()


def test_get_info_missing_argument():
    runner = CliRunner()

    # Assert the command raises a PipelineApi exception
    result = runner.invoke(pipelines_commands.pipelines, ["get-info"])

    # Assert the command failed due to missing argument
    assert result.exit_code != 0
    assert "Error: Missing argument 'PIPELINE_NAME'" in result.output


def test_get_info_api_exception(caplog, unstub):
    runner = CliRunner()

    when(pipelines_commands.pipelines_logic).get_pipeline_info(
        "bad_pipeline_name"
    ).thenRaise(
        ApiException(
            status=400,
            reason="Error Reason",
            body='{"message": "this is the body message"}',
        )
    )

    with caplog.at_level(logging.DEBUG):
        result = runner.invoke(
            pipelines_commands.pipelines, ["get-info", "bad_pipeline_name"]
        )

    # Assert the command failed and that the error handler formatted the error message
    assert result.exit_code != 0
    verify(pipelines_commands.pipelines_logic).get_pipeline_info("bad_pipeline_name")
    assert (
        "API call failed with status code 400 (Error Reason): this is the body message"
        in caplog.text
    )

    unstub()
