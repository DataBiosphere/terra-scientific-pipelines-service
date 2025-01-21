import json

from ..teaspoons_smoke_test_case import TeaspoonsSmokeTestCase


class TeaspoonsPipelinesListTests(TeaspoonsSmokeTestCase):
    @staticmethod
    def pipelines_list_url() -> str:
        return TeaspoonsSmokeTestCase.build_teaspoons_url("/api/pipelines/v1")

    def test_status_code_is_200(self):
        response = TeaspoonsSmokeTestCase.call_teaspoons(self.pipelines_list_url(), TeaspoonsSmokeTestCase.USER_TOKEN)
        self.assertEqual(response.status_code, 200, f"Pipelines List HTTP Status is not 200: {response.text}")

    def test_pipelines_list_greater_than_zero(self):
        response = TeaspoonsSmokeTestCase.call_teaspoons(self.pipelines_list_url(), TeaspoonsSmokeTestCase.USER_TOKEN)
        pipeline_info = json.loads(response.text)
        self.assertGreater(len(pipeline_info["results"]), 0, "No pipelines found")
