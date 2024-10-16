# utils.py

import click
import json

from pydantic import BaseModel


def _pretty_print(obj: BaseModel):
    """
    Prints a pydantic model in a pretty format to the console
    """
    click.echo(json.dumps(obj.model_dump(), indent=4))
