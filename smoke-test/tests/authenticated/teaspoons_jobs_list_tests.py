from ..teaspoons_smoke_test_case import TeaspoonsSmokeTestCase


class TeaspoonsJobsListTests(TeaspoonsSmokeTestCase):
    '''
    Test the jobs list endpoint for a 200 status code
    '''
    @staticmethod
    def jobs_list_url() -> str:
        return TeaspoonsSmokeTestCase.build_teaspoons_url("api/job/v1/jobs?limit=10")

    def test_status_code_is_200(self):
        response = TeaspoonsSmokeTestCase.call_teaspoons(self.jobs_list_url(), TeaspoonsSmokeTestCase.USER_TOKEN)
        self.assertEqual(response.status_code, 200, f"Jobs List HTTP Status is not 200: {response.text}")
