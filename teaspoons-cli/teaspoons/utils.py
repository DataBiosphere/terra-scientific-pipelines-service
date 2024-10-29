# utils.py

import json
import logging

from functools import wraps
from pydantic import BaseModel

from teaspoons_client.exceptions import ApiException


LOGGER = logging.getLogger(__name__)


def _pretty_print(obj: BaseModel, logger=LOGGER):
    """
    Prints a pydantic model in a pretty format to the console
    """
    try:
        LOGGER.info(json.dumps(obj.model_dump(), indent=4))
    except Exception:
        LOGGER.error(obj)


def handle_api_exceptions(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ApiException as e:
            formatted_message = f"API call failed with status code {e.status} ({e.reason}): {json.loads(e.body)['message']}"
            LOGGER.info(formatted_message)
            exit(1)
        except Exception as e:
            LOGGER.info(str(e))
            exit(1)

    return wrapper
