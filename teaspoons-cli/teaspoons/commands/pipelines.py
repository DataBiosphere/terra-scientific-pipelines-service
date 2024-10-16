import click
import typer

# import generated library
from generated.teaspoons_client import PipelinesApi

# teaspoons modules
from client import ClientWrapper
from utils import _pretty_print

pipelines_app = typer.Typer()


@pipelines_app.command()
def list():
    try:
        with ClientWrapper() as api_client:
            # NOTE: this is broken until this is merged with the different api def for the get_pipelines_result
            pipeline_client = PipelinesApi(api_client=api_client)
            pipelines = pipeline_client.get_pipelines()
            for pipeline in pipelines.pipelines.results:
                click.echo(f'{pipeline.name} - {pipeline.description}')
    except ValueError as e:
        click.echo(str(e), err=True)
        exit(1)


@pipelines_app.command()
@click.argument('name')
def get_info(name: str):
    try:
        with ClientWrapper() as api_client:
            pipeline_client = PipelinesApi(api_client=api_client)
            pipeline = pipeline_client.get_pipeline_details(pipeline_name=name)
            _pretty_print(pipeline)
    except ValueError as e:
        click.echo(str(e), err=True)
        exit(1)
