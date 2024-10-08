# config.py

from pathlib import Path
from oauth2_cli_auth import OAuth2ClientInfo
from dotenv import dotenv_values


def init():
    global config
    global client_info
    global token_file
    # setup
    home = Path.home()
    config = dotenv_values(".env")
    client_info = OAuth2ClientInfo.from_oidc_endpoint(
        config["OAUTH_OPENID_CONFIGURATION_URI"],
        client_id=config["OAUTH_CLIENT_ID"],
        scopes=[f"openid+email+profile+{config['OAUTH_CLIENT_ID']}"]
    )

    # Figure out the path to the token file...there must be some way to abstract out local storage
    if config["LOCAL_STORAGE_PATH"].startswith("/"):
        token_file = f"{config['LOCAL_STORAGE_PATH']}/access_token"
    else:
        token_file = f'{home}/{config["LOCAL_STORAGE_PATH"]}/access_token'
