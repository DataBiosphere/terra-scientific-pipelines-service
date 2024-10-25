# tests/logic/test_pipelines_logic.py

from unittest.mock import patch, Mock
import pytest
from teaspoons.logic import pipelines_logic
from teaspoons_client.exceptions import ApiException


@pytest.fixture
def mock_cli_config():
    with patch("teaspoons.logic.pipelines_logic.CliConfig") as mock_config:
        mock_config.token_file = "mock_token_file"
        yield mock_config


@pytest.fixture
def mock_client_wrapper():
    with patch("teaspoons.logic.pipelines_logic.ClientWrapper") as mock_wrapper:
        mock_client = Mock()
        mock_wrapper.return_value.__enter__.return_value = mock_client
        yield mock_client


@pytest.fixture
def mock_pipelines_api(mock_client_wrapper):
    with patch("teaspoons.logic.pipelines_logic.PipelinesApi") as mock_api:
        mock_instance = mock_api.return_value
        yield mock_instance


def test_list_pipelines(mock_pipelines_api):
    # Arrange
    mock_pipeline = Mock(pipeline_name="Test Pipeline", description="Test Description")
    mock_pipelines_api.get_pipelines.return_value.results = [mock_pipeline]

    # Act
    result = pipelines_logic.list_pipelines()

    # Assert
    assert len(result) == 1
    assert result[0].pipeline_name == "Test Pipeline"
    assert result[0].description == "Test Description"


def test_get_pipeline_info(mock_pipelines_api):
    # Arrange
    pipeline_name = "Test Pipeline"
    mock_pipeline_with_details = Mock(pipeline_name=pipeline_name, description="Test Description", display_name="Test Display Name", type="Test Type", inputs=[])
    mock_pipelines_api.get_pipeline_details.return_value = mock_pipeline_with_details

    # Act
    result = pipelines_logic.get_pipeline_info(pipeline_name)

    # Assert
    assert result == mock_pipeline_with_details


def test_get_pipeline_info_bad_pipeline_name(mock_pipelines_api):
    # Arrange
    pipeline_name = "Bad Pipeline Name"
    mock_pipelines_api.get_pipeline_details.side_effect = ApiException(
        404, reason="Pipeline not found"
    )

    # ApiException is raised
    with pytest.raises(ApiException):
        pipelines_logic.get_pipeline_info(pipeline_name)

    # Assert get_pipeline_details was called once
    mock_pipelines_api.get_pipeline_details.assert_called_once_with(pipeline_name=pipeline_name)
