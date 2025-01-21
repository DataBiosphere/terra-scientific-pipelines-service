import re
from functools import cache
from unittest import TestCase
from urllib.parse import urljoin

import requests
from requests import Response


class TeaspoonsSmokeTestCase(TestCase):
    TEASPOONS_HOST = None
    USER_TOKEN = None

    @staticmethod
    def build_teaspoons_url(path: str) -> str:
        assert TeaspoonsSmokeTestCase.TEASPOONS_HOST, "ERROR - TeaspoonsSmokeTests.TEASPOONS_HOST not properly set"
        if re.match(r"^\s*https?://", TeaspoonsSmokeTestCase.TEASPOONS_HOST):
            return urljoin(TeaspoonsSmokeTestCase.TEASPOONS_HOST, path)
        else:
            return urljoin(f"https://{TeaspoonsSmokeTestCase.TEASPOONS_HOST}", path)

    @staticmethod
    @cache
    def call_teaspoons(url: str, user_token: str = None) -> Response:
        """Function is memoized so that we only make the call once"""
        headers = {"Authorization": f"Bearer {user_token}"} if user_token else {}
        return requests.get(url, headers=headers)
