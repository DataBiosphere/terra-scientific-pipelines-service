# logic/pipelines_logic.py

from teaspoons_client import PipelinesApi

from client import ClientWrapper
from utils import _pretty_print, handle_api_exceptions


@handle_api_exceptions
def list_pipelines() -> None:
    with ClientWrapper() as api_client:
        pipeline_client = PipelinesApi(api_client=api_client)
        pipelines = pipeline_client.get_pipelines()

        for pipeline in pipelines.results:
            print(f"{pipeline.pipeline_name} - {pipeline.description}")


@handle_api_exceptions
def get_pipeline_info(pipeline_name: str) -> None:
    with ClientWrapper() as api_client:
        pipeline_client = PipelinesApi(api_client=api_client)
        pipeline = pipeline_client.get_pipeline_details(pipeline_name=pipeline_name)
        _pretty_print(pipeline)
