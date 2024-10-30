# cli.py

import click
import logging

from teaspoons import __version__, log
from teaspoons.commands.auth_commands import auth
from teaspoons.commands.pipelines_commands import pipelines


# Context settings for commands, for overwriting some click defaults
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

LOGGER = logging.getLogger(__name__)


@click.group(name="teaspoons", context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
@click.option(
    "--debug",
    is_flag=True,
    help="DEBUG-level logging",
)
def cli(debug):
    log.configure_logging(debug)
    LOGGER.debug(
        "Log level set to: %s", logging.getLevelName(logging.getLogger().level)
    )


cli.add_command(auth)
cli.add_command(pipelines)
# will add runs_app later


if __name__ == "__main__":
    cli()
