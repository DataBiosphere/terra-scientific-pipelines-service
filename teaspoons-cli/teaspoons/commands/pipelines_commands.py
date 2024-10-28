# commands/pipelines_commands.py

import click

from logic import pipelines_logic
from utils import _pretty_print, handle_api_exceptions


@click.group()
def pipelines():
    """Commands for running Teaspooons pipelines"""


@pipelines.command()
@handle_api_exceptions
def list():
    """List all available pipelines"""
    pipelines_list = pipelines_logic.list_pipelines()
    click.echo(f"Found {len(pipelines_list)} available pipeline{'' if len(pipelines_list) == 1 else 's'}:")
    for pipeline in pipelines_list:
        click.echo(f"\t{pipeline.pipeline_name} - {pipeline.description}")


@pipelines.command()
@click.argument("pipeline_name")
@handle_api_exceptions
def get_info(pipeline_name: str):
    """Get information about a specific pipeline"""
    pipeline_info = pipelines_logic.get_pipeline_info(pipeline_name)
    click.echo(_pretty_print(pipeline_info))
