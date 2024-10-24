# client.py

from generated.teaspoons_client import Configuration, ApiClient
from config import CliConfig
from auth_helper import _load_local_token


def _get_api_client(token: str, api_url: str) -> ApiClient:
    api_config = Configuration()
    api_config.host = api_url
    api_config.access_token = token
    return ApiClient(configuration=api_config)

class ClientWrapper:
    """
    Wrapper to ensure that the user is authenticated before running the callback and that provides the low level api client to be used
    by subsequent commands
    """

    def __enter__(self):
        cli_config = CliConfig()  # initialize the config from environment variables
        token = _load_local_token(cli_config.token_file)
        if not token:
            raise ValueError('Please authenticate first')
        else:
            return _get_api_client(token, cli_config.config["TEASPOONS_API_URL"])

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
