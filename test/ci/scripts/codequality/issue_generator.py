#!/usr/bin/env python3

from . import print_debug
import hashlib
import os


class IssueGenerator:
    """Interface for the issue generators"""

    def __init__(self, issue_list, **kwargs):
        self._issues = issue_list
        self._files = []

    def _format_issue(self, unfmt_issue):
        filepath = os.path.relpath(unfmt_issue['file'])
        issue = {
            'type': 'issue',
            'check_name': unfmt_issue['name'],
            'description': unfmt_issue['description'],
            'location': {
                "path": filepath,
                "lines": {
                    "begin": unfmt_issue['line'],
                    "end": unfmt_issue['line'],
                },
                "positions": {
                    "begin": {
                        "line": unfmt_issue['line'],
                        "column": unfmt_issue['column'],
                    },
                    'end': {
                        "line": unfmt_issue['line'],
                        "column": unfmt_issue['column'],
                    },
                },
            },
        }


        if 'end_line' in unfmt_issue:
            issue['location']['positions']['end'] = {
                "line": unfmt_issue['end_line'],
                "column": unfmt_issue['column'],
            }
            issue['location']['lines']['end'] = unfmt_issue['end_line']

        issue['fingerprint'] = hashlib.md5(
            '{file}:{line}:{column}:{type}'.format(
                file=filepath,
                line=unfmt_issue['line'],
                column=unfmt_issue['column'],
                type=unfmt_issue['name']).encode()).hexdigest()

        issue['categories'], issue['severity'] = \
            self._get_classifiaction(unfmt_issue)

        print_debug(issue)
        return issue
