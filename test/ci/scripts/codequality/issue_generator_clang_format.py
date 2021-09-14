#!/usr/bin/env python3

from . import print_debug, print_info
from .issue_generator_clang_tool import ClangToolIssueGenerator
import os
import re
import copy
import difflib
import subprocess


class ClangFormatIssueGenerator(ClangToolIssueGenerator):
    """issue generator for clang format"""

    def __init__(self, **kwargs):
        kwargs['clang_tool_executable'] = kwargs.pop('clang_format_executable',
                                                     'clang-format')
        super().__init__('clang-format', **kwargs)

    def _get_classifiaction(self, issue):
        return (['Style'], 'info')

    def generate_issues(self):
        issue = {}
        for filename in self._files:
            with open(filename, 'r') as fh:
                unformated_file = fh.readlines()

            command = copy.copy(self._command)
            command.append(filename)
            formated_file = list(self._run_command(command))

            # diffs = difflib.unified_diff(unformated_file, formated_file, n=0)

            # print(diffs)
            # for diff in diffs:
            #     print(diff, end='')

            s = difflib.SequenceMatcher(None, unformated_file, formated_file)
            for tag, i1, i2, j1, j2 in s.get_opcodes():
                if tag != 'equal':
                    diff = list(
                        difflib.unified_diff(
                            unformated_file[i1:i2],
                            formated_file[j1:j2]))
                    issue = {
                        'name': f'''clang-format:{tag}''',
                        'description': ''.join(diff[3:]),
                        'file': filename,
                        'line': i1,
                        'column': 1,
                        'end_line': i2,
                    }
                    self.add_issue(issue)
