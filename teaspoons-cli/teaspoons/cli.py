# cli.py

import typer
from commands.auth import auth_app
from commands.pipelines import pipelines_app

cli = typer.Typer()
cli.add_typer(auth_app, name="auth", help="Authentication commands")
cli.add_typer(pipelines_app, name="pipelines", help="Pipeline commands")
# will add runs_app later


if __name__ == '__main__':
    cli()
