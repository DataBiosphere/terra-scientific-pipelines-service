# commands/auth_commands.py

import click

from logic import auth_logic

@click.group()
def auth():
    """Commands for authenticating to Teaspoons"""


@auth.command()
def login():
    """Authenticate with Teaspoons via browser login to Terra b2c"""
    auth_logic.check_local_token_and_fetch_if_needed()


@auth.command()
def logout():
    """Clear the local authentication token"""
    auth_logic.clear_local_token()
