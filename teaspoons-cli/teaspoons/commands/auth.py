# auth.py

import typer
import click

from auth_helper import get_access_token_with_browser_open, _validate_token, _save_local_token, _load_local_token, _clear_local_token
from config import CliConfig

cli_config = CliConfig()  # initialize the config from environment variables

auth_app = typer.Typer()


@auth_app.command()
def login():
    token = _load_local_token(cli_config.token_file)
    if token and _validate_token(token):
        click.echo('Already authenticated')
        return
    token = get_access_token_with_browser_open(cli_config.client_info)
    _save_local_token(cli_config.token_file, token)


@auth_app.command()
def logout():
    _clear_local_token(cli_config.token_file)
    click.echo('Logged out')
