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
import subprocess
try:
    from termcolor import colored
except ImportError:
    def colored(text, color):
        return text


def print_debug(message):
    '''helper finction to print debug messages'''
    print(f'Debug: {colored(message, "red")}',
          file=sys.stderr, flush=True)


class ClangTidy2CodeQuality:
    '''
    Main class to run and convert the results of clang-tidy to the code-quality
    format
    '''
    # 7-bit C1 ANSI sequences
    ANSI_ESCAPE = re.compile(r'''
        \x1B  # ESC
        (?:   # 7-bit C1 Fe (except CSI)
            [@-Z\\-_]
        |     # or [ for CSI, followed by a control sequence
            \[
            [0-?]*  # Parameter bytes
            [ -/]*  # Intermediate bytes
            [@-~]   # Final byte
        )
    ''', re.VERBOSE)

    ISSUE_PARSE = re.compile(r'(?P<file>.*\.(cc|hh)):(?P<line>[0-9]+):(?P<column>[0-9]+): (warning|error): (?P<detail>.*) \[(?P<type>.*)\]')  # noqa

    CLASSIFICATIONS = {
        'bugprone': {
            'categories': ['Bug Risk'],
            'severity': 'major',
        },
        'modernize': {
            'categories': ['Clarity', 'Compatibility', 'Style'],
            'severity': 'info'
        },
        'mpi': {
            'categories': ['Bug Risk', 'Performance'],
            'severity': 'critical',
        },
        'openmp': {
            'categories': ['Bug Risk', 'Performance'],
            'severity': 'critical',
        },
        'performance': {
            'categories': ['Performance'],
            'severity': 'minor',
        },
        'readability': {
            'categories': ['Clarity', 'Style'],
            'severity': 'info'
        },
    }

    def __init__(self):
        self._issues = {}
        self._command = ['run-clang-tidy']
        self._command.extend(sys.argv[1:])

    def run(self):
        '''run clang tidy and generage a code quality report'''
        print_debug(f'[clang-tidy] command: {self._command}')
        self._generate_issues(self._command)
        print(json.dumps(list(self._issues.values())))

    @property
    def command(self):
        '''get the command that is run'''
        return self._command

    def _get_classifiaction(self, type_):
        categories = ['Bug Risk']
        severity = 'blocker'

        if type_ in self.CLASSIFICATIONS:
            categories = self.CLASSIFICATIONS[type_]['categories']
            severity = self.CLASSIFICATIONS[type_]['severity']
        elif type_[0] == 'clang':
            if type_[1] == 'diagnostic':
                categories = ['Bug Risk']
                severity = 'blocker'
            elif type_[1] == 'analyzer':
                categories = ['Bug Risk']
                severity = 'major'

        return (categories, severity)

    def _run_command(self, command):
        popen = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 universal_newlines=True)

        for stdout_line in iter(popen.stdout.readline, ""):
            clean_line = self.ANSI_ESCAPE.sub('', stdout_line).rstrip()
            print_debug(clean_line)
            yield clean_line

        popen.stdout.close()

        return_code = popen.wait()
        if return_code:
            print_debug(
                f"[clang-tidy] {command} ReturnCode {return_code}")

    def _generate_issues(self, command):
        issue = {}
        for line in self._run_command(command):
            match = self.ISSUE_PARSE.match(line)
            if match:
                if len(issue) != 0:
                    self._add_issue(issue)
                issue = match.groupdict()
            elif issue:
                if 'content' in issue:
                    issue['content'].append(line)
                else:
                    issue['content'] = [line]
        self._add_issue(issue)

    def _add_issue(self, issue):
        issue_ = self._format_issue(issue)

        if issue_['fingerprint'] in self._issues:
            return

        self._issues[issue_['fingerprint']] = issue_

    def _format_issue(self, issue_dict):
        issue_dict['file'] = os.path.relpath(issue_dict['file'])

        issue = {
            'type': 'issue',
            'check_name': issue_dict['type'],
            'description': issue_dict['detail'],
            'location': {
                "path": issue_dict['file'],
                "lines": {
                    "begin": int(issue_dict['line']),
                    "end": int(issue_dict['line']),
                },
                "positions": {
                    "begin": {
                        "line": int(issue_dict['line']),
                        "column": int(issue_dict['column']),
                    },
                    "end": {
                        "line": int(issue_dict['line']),
                        "column": int(issue_dict['column']),
                    },
                },
            },
        }

        if 'content' in issue_dict:
            issue['content'] = {
                'body': '```\n' +
                '\n'.join(issue_dict['content']) +
                '\n```'
            }

        issue['fingerprint'] = hashlib.md5(
            '{file}:{line}:{column}:{type}'.format(**issue_dict).encode()
        ).hexdigest()

        type_ = issue_dict['type'].split('-')[0]
        issue['categories'], issue['severity'] = self._get_classifiaction(type_)

        return issue


formater = ClangTidy2CodeQuality()
formater.run()
