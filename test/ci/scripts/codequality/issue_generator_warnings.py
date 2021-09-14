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

from . import print_debug, print_info
from .issue_generator import IssueGenerator
import re
import sys
import warning_parser as warn


class WarningsIssueGenerator(IssueGenerator):
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        files = kwargs.pop('files')

        self._input_files = {}
        compiler_re = re.compile(".*build.*(gcc|clang)-err\.log")

        for _file in files:
            match = compiler_re.search(_file)
            if match:
                self._input_files[match.group(1)] = _file
            else:
                print_info(f"Skipped {_file}, could not determine compiler")

    def generate_issues(self):
        '''parse warning files'''

        for compiler, _file in self._input_files.items():
            warnings = warn.get_warnings(_file, compiler)
            for warning in warnings:
                issue = {
                    'name': f'warning:{compiler}-{warning.get_category()}',
                    'description': warning.get_message(),
                    'file': warning.get_filepath(),
                    'line': warning.get_line(),
                    'column': warning.get_column(),
                    'raw': warning,
                }

                self.add_issue(issue)

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
