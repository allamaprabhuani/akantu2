#!/usr/bin/env python3
"""clang-tidy2code-quality.py: Conversion of clang-tidy output 2
code-quality"""

__author__ = "Nicolas Richart"
__credits__ = [
    "Nicolas Richart <nicolas.richart@epfl.ch>",
]
__copyright__ = "Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"

import re
import os
import hashlib
import json
import sys
import warning_parser as warn
try:
    from termcolor import colored
except ImportError:
    def colored(text, color):
        '''replace colored if not present'''
        return text


def print_debug(message):
    '''helper finction to print debug messages'''
    print(f'Debug: {colored(message, "red")}',
          file=sys.stderr, flush=True)


class Warnings2CodeQuality:
    '''
    Main class to run and convert the results of clang-tidy to the code-quality
    format
    '''
    CLASSIFICATIONS = {
        'gcc': {
            'uninitialized': {
                'categories': ['Bug Risk'],
                'severity': 'major',
            },
            'sign-compare': {
                'categories': ['Bug Risk'],
                'severity': 'minor'
            },
        },
    }

    def __init__(self):
        self._issues = {}
        self._files = sys.argv[1:]

    def parse(self):
        '''parse warning files'''
        compiler_re = re.compile(".*build.*(gcc|clang)-err\.log")

        for _file in self._files:
            match = compiler_re.search(_file)
            if not match:
                continue

            compiler = match.group(1)
            warnings = warn.get_warnings(_file, compiler)
            self._add_issues(warnings)

        print(json.dumps(list(self._issues.values())))

    @property
    def list_files(self):
        '''get the list of files to analyse'''
        return self._files

    def _get_classifiaction(self, warning):
        categories = ['Clarity']
        severity = 'info'

        if warning.get_tool() in self.CLASSIFICATIONS:
            classifications = self.CLASSIFICATIONS[warning.get_tool()]
            if warning.get_category() in classifications:
                cat = warning.get_category()
                categories = classifications[cat]['categories']
                severity = classifications[cat]['severity']

        return (categories, severity)

    def _add_issues(self, warnings):
        for warning in warnings:
            issue_ = self._format_issue(warning)

            if issue_['fingerprint'] in self._issues:
                continue

            print_debug(f'{issue_}')
            self._issues[issue_['fingerprint']] = issue_

    def _format_issue(self, warning):
        issue = {
            'type': 'issue',
            'check_name': warning.get_category(),
            'description': warning.get_message(),
            'location': {
                "path": warning.get_filepath(),
                "lines": {
                    "begin": warning.get_line(),
                    "end": warning.get_line(),
                },
                "positions": {
                    "begin": {
                        "line": warning.get_line(),
                        "column": warning.get_column(),
                    },
                    "end": {
                        "line": warning.get_line(),
                        "column": warning.get_column(),
                    },
                },
            },
        }

        issue['fingerprint'] = hashlib.md5(
            '{file}:{line}:{column}:{type}'.format(
                file=warning.get_filepath(),
                line=warning.get_line(),
                column=warning.get_column(),
                type=warning.get_category()).encode()
        ).hexdigest()

        issue['categories'], issue['severity'] = self._get_classifiaction(warning)

        return issue


formater = Warnings2CodeQuality()
formater.parse()
