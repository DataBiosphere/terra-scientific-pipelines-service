# tests/logic/test_pipelines_logic.py

from unittest.mock import patch, Mock
import pytest
from teaspoons.logic import pipelines_logic
from teaspoons_client.exceptions import ApiException


@pytest.fixture
def mock_cli_config():
    with patch('teaspoons.logic.pipelines_logic.CliConfig') as mock_config:
        mock_config.token_file = 'mock_token_file'
        yield mock_config

@pytest.fixture
def mock_client_wrapper():
    with patch('teaspoons.logic.pipelines_logic.ClientWrapper') as mock_wrapper:
        mock_client = Mock()
        mock_wrapper.return_value.__enter__.return_value = mock_client
        yield mock_client

@pytest.fixture
def mock_pipelines_api(mock_client_wrapper):
    with patch('teaspoons.logic.pipelines_logic.PipelinesApi') as mock_api:
        mock_instance = mock_api.return_value
        yield mock_instance


def test_list_pipelines(mock_pipelines_api, capsys):
    # Arrange
    mock_pipeline = Mock(pipeline_name="Test Pipeline", description="Test Description")
    mock_pipelines_api.get_pipelines.return_value.results = [mock_pipeline]

    # Act
    pipelines_logic.list_pipelines()

    # Assert
    captured = capsys.readouterr()
    mock_pipelines_api.get_pipelines.assert_called_once()
    assert captured.out == "Test Pipeline - Test Description\n"


def test_get_pipeline_info(mock_pipelines_api, capsys):
    # Arrange
    pipeline_name = "Test Pipeline"
    expected_output = {"pipeline_name": pipeline_name, "description": "Test Description"}
    mock_pipelines_api.get_pipeline_details.return_value = expected_output

    # Act
    pipelines_logic.get_pipeline_info(pipeline_name)

    # Assert
    captured = capsys.readouterr()
    mock_pipelines_api.get_pipeline_details.assert_called_once_with(pipeline_name=pipeline_name)
    assert captured.out == f"{expected_output}\n"


def test_get_pipeline_info_bad_pipeline_name(mock_pipelines_api, capsys):
    # Arrange
    pipeline_name = "Bad Pipeline Name"
    mock_pipelines_api.get_pipeline_details.side_effect = ApiException(404, reason="Not Found", body="{\"message\": \"Pipeline not found\"}")

    # Act
    with pytest.raises(SystemExit) as e:
        pipelines_logic.get_pipeline_info(pipeline_name)

    # Assert
    captured = capsys.readouterr()
    assert captured.out == "API call failed with status code 404 (Not Found): Pipeline not found\n"
