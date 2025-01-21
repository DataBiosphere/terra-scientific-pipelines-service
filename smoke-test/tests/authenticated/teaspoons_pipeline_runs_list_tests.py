from ..teaspoons_smoke_test_case import TeaspoonsSmokeTestCase


class TeaspoonsPipelineRunsListTests(TeaspoonsSmokeTestCase):
    @staticmethod
    def pipeline_runs_list_url() -> str:
        return TeaspoonsSmokeTestCase.build_teaspoons_url("api/pipelineruns/v1/pipelineruns?limit=10")

    def test_status_code_is_200(self):
        response = TeaspoonsSmokeTestCase.call_teaspoons(self.pipeline_runs_list_url(), TeaspoonsSmokeTestCase.USER_TOKEN)
        self.assertEqual(response.status_code, 200, f"Pipeline Runs List HTTP Status is not 200: {response.text}")
