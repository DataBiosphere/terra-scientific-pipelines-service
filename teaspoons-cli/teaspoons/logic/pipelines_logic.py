# logic/pipelines_logic.py

from teaspoons_client import PipelinesApi
from teaspoons_client.models.pipeline_with_details import PipelineWithDetails
from teaspoons_client.models.pipeline import Pipeline

from client import ClientWrapper


def list_pipelines() -> list[Pipeline]:
    """List all pipelines, returning a list of dictionaries."""
    with ClientWrapper() as api_client:
        pipeline_client = PipelinesApi(api_client=api_client)
        pipelines = pipeline_client.get_pipelines()

        result = []
        for pipeline in pipelines.results:
            result.append(pipeline)

        return result


def get_pipeline_info(pipeline_name: str) -> PipelineWithDetails:
    """Get the details of a pipeline, returning a dictionary."""
    with ClientWrapper() as api_client:
        pipeline_client = PipelinesApi(api_client=api_client)
        return pipeline_client.get_pipeline_details(pipeline_name=pipeline_name)
