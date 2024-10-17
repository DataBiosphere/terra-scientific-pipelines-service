# config.py

from pathlib import Path
from oauth2_cli_auth import OAuth2ClientInfo
from dotenv import dotenv_values


class CliConfig:
    """A class to hold configuration information for the CLI"""

    def __init__(self):
        self.config = dotenv_values("../.env")

        self.client_info = OAuth2ClientInfo.from_oidc_endpoint(
            self.config["OAUTH_OPENID_CONFIGURATION_URI"],
            client_id=self.config["OAUTH_CLIENT_ID"],
            scopes=[f"openid+email+profile+{self.config['OAUTH_CLIENT_ID']}"]
        )

        self.server_port = int(self.config["SERVER_PORT"])

        self.token_file = f'{Path.home()}/{self.config["LOCAL_STORAGE_PATH"]}/access_token'
