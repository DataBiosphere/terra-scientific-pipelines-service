# pipelines.py

import click
import typer

# import generated library
from teaspoons_client import PipelinesApi
from utils import handle_exceptions

# teaspoons modules
from client import ClientWrapper
from utils import _pretty_print

pipelines_app = typer.Typer()


@pipelines_app.command()
@handle_exceptions
def list():
    with ClientWrapper() as api_client:
        pipeline_client = PipelinesApi(api_client=api_client)
        pipelines = pipeline_client.get_pipelines()
        for pipeline in pipelines.pipelines.results:
            click.echo(f'{pipeline.name} - {pipeline.description}')

@pipelines_app.command()
@click.argument('name')
@handle_exceptions
def get_info(name: str):
    with ClientWrapper() as api_client:
        pipeline_client = PipelinesApi(api_client=api_client)
        pipeline = pipeline_client.get_pipeline_details(pipeline_name=name)
        _pretty_print(pipeline)
