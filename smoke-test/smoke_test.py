import argparse
import sys
import unittest
from unittest import TestSuite

import requests

from tests.authenticated.teaspoons_jobs_list_tests import TeaspoonsJobsListTests
from tests.authenticated.teaspoons_pipeline_runs_list_tests import TeaspoonsPipelineRunsListTests
from tests.authenticated.teaspoons_pipelines_list_tests import TeaspoonsPipelinesListTests
from tests.teaspoons_smoke_test_case import TeaspoonsSmokeTestCase
from tests.unauthenticated.teaspoons_status_tests import TeaspoonsStatusTests
from tests.unauthenticated.teaspoons_version_tests import TeaspoonsVersionTests

DESCRIPTION = """
Teaspoons Smoke Test
Enter the host (domain and optional port) of the Teaspoons instance you want to to test. 
This test will ensure that the Teaspoons instance running on that host is minimally functional.
"""


def gather_tests(is_authenticated: bool = False) -> TestSuite:
    suite = unittest.TestSuite()

    status_tests = unittest.defaultTestLoader.loadTestsFromTestCase(TeaspoonsStatusTests)
    version_tests = unittest.defaultTestLoader.loadTestsFromTestCase(TeaspoonsVersionTests)

    suite.addTests(status_tests)
    suite.addTests(version_tests)

    if is_authenticated:
        pipeline_list_tests = unittest.defaultTestLoader.loadTestsFromTestCase(TeaspoonsPipelinesListTests)
        pipeline_runs_list_tests = unittest.defaultTestLoader.loadTestsFromTestCase(TeaspoonsPipelineRunsListTests)
        jobs_list_tests = unittest.defaultTestLoader.loadTestsFromTestCase(TeaspoonsJobsListTests)

        suite.addTests(pipeline_list_tests)
        suite.addTests(pipeline_runs_list_tests)
        suite.addTests(jobs_list_tests)
    else:
        print("No User Token provided.  Skipping authenticated tests.")

    return suite


def main(main_args):
    if main_args.user_token:
        verify_user_token(main_args.user_token)

    TeaspoonsSmokeTestCase.TEASPOONS_HOST = main_args.teaspoons_host
    TeaspoonsSmokeTestCase.USER_TOKEN = main_args.user_token

    test_suite = gather_tests(main_args.user_token)

    runner = unittest.TextTestRunner(verbosity=main_args.verbosity)
    result = runner.run(test_suite)

    # system exit if any tests fail
    if result.failures or result.errors:
        sys.exit(1)


def verify_user_token(user_token: str) -> bool:
    response = requests.get(f"https://www.googleapis.com/oauth2/v1/tokeninfo?access_token={user_token}")
    assert response.status_code == 200, "User Token is no longer valid.  Please generate a new token and try again."


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(
            description=DESCRIPTION,
            formatter_class=argparse.RawTextHelpFormatter
        )
        parser.add_argument(
            "-v",
            "--verbosity",
            type=int,
            choices=[0, 1, 2],
            default=1,
            help="""Python unittest verbosity setting: 
0: Quiet - Prints only number of tests executed
1: Minimal - (default) Prints number of tests executed plus a dot for each success and an F for each failure
2: Verbose - Help string and its result will be printed for each test"""
        )
        parser.add_argument(
            "teaspoons_host",
            help="domain with optional port number of the Teasponns host you want to test"
        )
        parser.add_argument(
            "user_token",
            nargs='?',
            default=None,
            help="Optional. If present, will test additional authenticated endpoints using the specified token"
        )

        args = parser.parse_args()

        # Need to pop off sys.argv values to avoid messing with args passed to unittest.main()
        for _ in range(len(sys.argv[1:])):
            sys.argv.pop()

        main(args)
        sys.exit(0)

    except Exception as e:
        print(e)
        sys.exit(1)
