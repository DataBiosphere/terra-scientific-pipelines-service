
import json
import click
from pydantic import BaseModel

# import generated library
from generated.teaspoons_client import Configuration, ApiClient, PipelinesApi

# teaspoons modules
import config
from auth import get_access_token_with_browser_open, __validate_token, __save_local_token, __load_local_token, __clear_local_token


cli_config = config.CliConfig()  # initialize the config from environment variables


def __get_api_client(token: str) -> ApiClient:
    api_config = Configuration()
    api_config.host = cli_config.config["TEASPOONS_API_URL"]
    api_config.access_token = token
    return ApiClient(configuration=api_config)


class ClientWrapper(object):
    """
    Wrapper to ensure that the user is authenticated before running the callback and that provides the low level api client to be used
    by subsequent commands
    """

    def __enter__(self):
        token = globals()['__load_local_token'](cli_config.token_file)
        if not token:
            click.echo('Please authenticate first')
            # Don't love exiting here, but it's the easiest way to ensure the user is prompted to authenticate without getting a whole
            # stack trace
            exit(1)
        else:
            return globals()['__get_api_client'](token)

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

def __pretty_print(obj: BaseModel):
    """
    Prints a pydantic model in a pretty format to the console
    """
    click.echo(json.dumps(obj.model_dump(), indent=4))


# CLI definition
@click.group()
def cli():
    pass


@cli.command()
@click.option('-t', '--token', default=False, is_flag=True, help='If specified, will ask you to provide a token manually')
def login(token: bool):
    if (token):
        # This is work around a weird behavior where python doesn't let you paste more than 1024 characters into the terminal.
        # Just importing fixes the issue...

        token = click.prompt('Please enter your access token', hide_input=False)
        if token and __validate_token(token):
            __save_local_token(cli_config.token_file, token)
            return
        else:
            click.echo('Invalid token')
            return

    else:
        token = __load_local_token(cli_config.token_file)
        if token:
            click.echo('Already authenticated')
            return
        token = get_access_token_with_browser_open(cli_config.client_info, server_port=10444)
        __save_local_token(cli_config.token_file, token)


@cli.command()
def logout():
    __clear_local_token(cli_config.token_file)
    click.echo('Logged out')


@cli.command()
def list_pipelines():
    with ClientWrapper() as api_client:
        # NOTE: this is broken until this is merged with the different api def for the get_pipelines_result
        pipeline_client = PipelinesApi(api_client=api_client)
        pipelines = pipeline_client.get_pipelines()
        for pipeline in pipelines.pipelines.results:
            click.echo(f'{pipeline.name} - {pipeline.description}')


@cli.command()
@click.argument('name')
def get_pipeline(name: str):
    with ClientWrapper() as api_client:
        pipeline_client = PipelinesApi(api_client=api_client)
        pipeline = pipeline_client.get_pipeline_details(pipeline_name=name)
        __pretty_print(pipeline)

