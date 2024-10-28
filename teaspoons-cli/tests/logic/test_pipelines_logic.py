# tests/logic/test_pipelines_logic.py

import pytest
from mockito import when, mock, verify

from teaspoons.logic import pipelines_logic
from teaspoons_client.exceptions import ApiException


@pytest.fixture
def mock_cli_config(unstub):
    config = mock({'token_file': 'mock_token_file'})
    when(pipelines_logic).CliConfig(...).thenReturn(config)
    yield config
    unstub()


@pytest.fixture
def mock_client_wrapper(unstub):
    client = mock()
    # Make the mock support context manager protocol
    when(client).__enter__().thenReturn(client)
    when(client).__exit__(None, None, None).thenReturn(None)
    
    when(pipelines_logic).ClientWrapper(...).thenReturn(client)
    yield client
    unstub()


@pytest.fixture
def mock_pipelines_api(mock_client_wrapper, unstub):
    api = mock()
    when(pipelines_logic).PipelinesApi(...).thenReturn(api)
    yield api
    unstub()


def test_list_pipelines(mock_pipelines_api, unstub):
    # Arrange
    mock_pipeline = mock({
        'pipeline_name': 'Test Pipeline',
        'description': 'Test Description'
    })
    when(mock_pipelines_api).get_pipelines().thenReturn(mock({'results': [mock_pipeline]}))

    # Act
    result = pipelines_logic.list_pipelines()

    # Assert
    assert len(result) == 1
    assert result[0].pipeline_name == "Test Pipeline"
    assert result[0].description == "Test Description"


def test_get_pipeline_info(mock_pipelines_api, unstub):
    # Arrange
    pipeline_name = "Test Pipeline"
    mock_pipeline = mock({
        'pipeline_name': pipeline_name,
        'description': 'Test Description',
        'display_name': 'Test Display Name',
        'type': 'Test Type',
        'inputs': []
    })
    when(mock_pipelines_api).get_pipeline_details(pipeline_name=pipeline_name).thenReturn(mock_pipeline)

    # Act
    result = pipelines_logic.get_pipeline_info(pipeline_name)

    # Assert
    assert result == mock_pipeline


def test_get_pipeline_info_bad_pipeline_name(mock_pipelines_api, unstub):
    # Arrange
    pipeline_name = "Bad Pipeline Name"
    when(mock_pipelines_api).get_pipeline_details(pipeline_name=pipeline_name).thenRaise(
        ApiException(404, reason="Pipeline not found")
    )

    # Act & Assert
    with pytest.raises(ApiException):
        pipelines_logic.get_pipeline_info(pipeline_name)

    # Verify the method was called
    verify(mock_pipelines_api).get_pipeline_details(pipeline_name=pipeline_name)
