# utils.py

import click
import json

from functools import wraps
from pydantic import BaseModel

from teaspoons_client.exceptions import ApiException


def _pretty_print(obj: BaseModel):
    """
    Prints a pydantic model in a pretty format to the console
    """
    click.echo(json.dumps(obj.model_dump(), indent=4))


def handle_exceptions(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ApiException as e:
            formatted_message = f"API call failed with status code {e.status} ({e.reason}): {json.loads(e.body)['message']}"
            click.echo(formatted_message, err=True)
            exit(1)
        except ValueError as e:
            click.echo(str(e), err=True)
            exit(1)
    return wrapper
