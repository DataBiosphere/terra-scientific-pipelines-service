# tests/test_config

from mockito import mock, when
from teaspoons import config


def test_config():
    mock_client_info = mock()
    when(config.OAuth2ClientInfo).from_oidc_endpoint(...).thenReturn(mock_client_info)

    test_config = config.CliConfig("tests/.test.config")

    assert test_config.config["TEASPOONS_API_URL"] == "not-real"
    assert test_config.config["SERVER_PORT"] == "12345"
    assert test_config.config["OAUTH_OPENID_CONFIGURATION_URI"] == "https://dontcare"
    assert test_config.config["OAUTH_CLIENT_ID"] == "whatever"
    assert test_config.config["LOCAL_STORAGE_PATH"] == ".cool"

    assert test_config.server_port == 12345
    assert test_config.client_info == mock_client_info
