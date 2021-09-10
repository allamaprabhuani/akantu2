#!/usr/bin/env python3
from . import print_debug, print_info
import os
import re


class IssueList:
    """Handles list of CodeQuality issues"""

    def __init__(self, extensions=None, excludes=None, **kwargs):
        if excludes is None:
            excludes = []

        if extensions is None:
            extensions = ['.cc', '.hh']

        self._extensions = [
            re.compile(r"\{}$".format(extension)) for extension in extensions
        ]

        self._exclude_patterns = [
            re.compile(exclude) for exclude in excludes
        ]

        self._issues = {}

    def _need_exclude(self, filename):
        need_exclude = False
        for pattern in self._exclude_patterns:
            match = pattern.search(filename)
            need_exclude |= bool(match)

        match_extension = False
        for extension in self._extensions:
            match = extension.search(filename)
            match_extension |= bool(match)

        need_exclude |= not match_extension

    def add_issue(self, issue):
        """add an issue to the list if not already present"""
        filepath = issue['location']['path']
        if self._need_exclude(filepath):
            return

        if issue['fingerprint'] in self._issues:
            return

        self._issues[issue['fingerprint']] = issue

    @property
    def issues(self):
        """get the list of registered issues"""
        return list(self._issues.values())
