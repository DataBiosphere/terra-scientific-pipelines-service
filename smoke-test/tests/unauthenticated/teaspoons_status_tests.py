from ..teaspoons_smoke_test_case import TeaspoonsSmokeTestCase


class TeaspoonsStatusTests(TeaspoonsSmokeTestCase):
    '''
    Test the status endpoint for a 200 status code
    '''
    @staticmethod
    def status_url() -> str:
        return TeaspoonsSmokeTestCase.build_teaspoons_url("/status")

    def test_status_code_is_200(self):
        response = TeaspoonsSmokeTestCase.call_teaspoons(self.status_url())
        self.assertEqual(response.status_code, 200)
