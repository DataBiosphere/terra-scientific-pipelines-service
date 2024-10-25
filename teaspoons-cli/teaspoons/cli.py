# cli.py

import click

from teaspoons import __version__
from teaspoons.commands.auth_commands import auth as auth
from teaspoons.commands.pipelines_commands import pipelines as pipelines


# Context settings for commands, for overwriting some click defaults
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="teaspoons", context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
def cli():
    pass


cli.add_command(auth)
cli.add_command(pipelines)
# will add runs_app later


if __name__ == "__main__":
    cli()
