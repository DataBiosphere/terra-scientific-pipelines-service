# tests/test_utils.py

import pytest
from unittest.mock import Mock, patch


@pytest.fixture
def mock_client_wrapper():
    with patch('teaspoons.client.ClientWrapper') as mock_client_wrapper:
        mock_client = Mock()
        mock_client_wrapper.return_value.__enter__.return_value = mock_client
        yield mock_client
