import json

from ..teaspoons_smoke_test_case import TeaspoonsSmokeTestCase


class TeaspoonsVersionTests(TeaspoonsSmokeTestCase):
    '''
    Test the version endpoint for a 200 status code and that 'build' is populated in the response
    '''
    @staticmethod
    def version_url() -> str:
        return TeaspoonsSmokeTestCase.build_teaspoons_url("/version")

    def test_status_code_is_200(self):
        response = TeaspoonsSmokeTestCase.call_teaspoons(self.version_url())
        self.assertEqual(response.status_code, 200)

    def test_version_value_specified(self):
        response = TeaspoonsSmokeTestCase.call_teaspoons(self.version_url())
        version = json.loads(response.text)
        self.assertIsNotNone(version["build"], "build value must be non-empty")
