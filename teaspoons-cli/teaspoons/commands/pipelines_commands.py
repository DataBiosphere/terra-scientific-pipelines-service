# commands/pipelines_commands.py

import click

from logic import pipelines_logic


@click.group()
def pipelines():
    """Commands for running Teaspooons pipelines"""


@pipelines.command()
def list():
    """List all available pipelines"""
    pipelines_logic.list_pipelines()


@pipelines.command()
@click.argument("pipeline_name")
def get_info(pipeline_name: str):
    """Get information about a specific pipeline"""
    pipelines_logic.get_pipeline_info(pipeline_name)
