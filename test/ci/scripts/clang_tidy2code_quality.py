#!/usr/bin/env python3
"""clang_tidy2code_quality.py: Conversion of clang-tidy output 2
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
import argparse
import copy
try:
    from termcolor import colored
except ImportError:
    def colored(text, color):  # pylint: disable=unused-argument
        """fallback function for termcolor.colored"""
        return text


def print_debug(message):
    '''helper function to print debug messages'''
    print(f'Debug: {colored(message, "red")}',
          file=sys.stderr, flush=True)


def print_info(message):
    '''helper function to print info messages'''
    print(f'Info: {colored(message, "blue")}',
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

    ISSUE_PARSE = re.compile(r'(?P<file>.*\.(cc|hh)):(?P<line>[0-9]+):(?P<column>[0-9]+): (warning|error): (?P<detail>.*) \[(?P<type>.*)\]')  # NOQA pylint: disable=line-too-long

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

    def __init__(self, compiledb_path, clang_tidy='clang-tidy', **kwargs):
        arguments = kwargs.pop('arguments', None)
        excludes = kwargs.pop('excludes', None)
        file_list = kwargs.pop('file_list', None)

        self._issues = {}
        self._command = [clang_tidy, '-p', compiledb_path]
        if arguments is None:
            arguments = []
        if excludes is None:
            excludes = []

        self._extensions = [re.compile(r'\.cc$'), re.compile(r'\.hh$')]

        self._command.extend(arguments)

        self._files = []

        self._exclude_patterns = []
        for exclude in excludes:
            self._exclude_patterns.append(re.compile(exclude))

        if file_list is None:
            file_list = self._get_files_from_compile_db(compiledb_path)

        self._define_file_list(file_list)

    def _get_files_from_compile_db(self, compiledb_path):  # pylint: disable=no-self-use
        file_list = []
        with open(os.path.join(
                compiledb_path,
                'compile_commands.json'), 'r') as compiledb_fh:
            compiledb = json.load(compiledb_fh)
            for entry in compiledb:
                file_list.append(entry['file'])
        return file_list

    def _define_file_list(self, file_list):
        for filename in file_list:
            filename = os.path.relpath(filename)
            need_exclude = self._need_exclude(filename)
            if need_exclude:
                print_debug(f'[clang-tidy] exluding file: {filename}')
                continue
            print_info(f'[clang-tidy] adding file: {filename}')
            self._files.append(filename)

    def run(self):
        '''run clang tidy and generage a code quality report'''
        self._generate_issues()
        print(json.dumps(list(self._issues.values())))

    @property
    def command(self):
        '''get the command that is run'''
        return self._command

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

        return need_exclude

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
        print_info(f'''[clang-tidy] command: {' '.join(command)}''')
        popen = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.DEVNULL,
                                 universal_newlines=True)

        for stdout_line in iter(popen.stdout.readline, ""):
            clean_line = self.ANSI_ESCAPE.sub('', stdout_line).rstrip()
            yield clean_line

        popen.stdout.close()

        return_code = popen.wait()
        if return_code:
            print_debug(
                f"[clang-tidy] {command} ReturnCode {return_code}")

    def _generate_issues(self):
        issue = {}
        for filename in self._files:
            command = copy.copy(self._command)
            command.append(filename)
            for line in self._run_command(command):
                match = self.ISSUE_PARSE.match(line)
                if match:
                    if len(issue) != 0:
                        self._add_issue(issue)
                    issue = match.groupdict()
                    print_debug(f'[clang-tidy] new issue: {line}')
                elif issue:
                    if 'content' in issue:
                        issue['content'].append(line)
                        print_debug(f'[clang-tidy] more extra content: {line}')
                    else:
                        issue['content'] = [line]
                        print_debug(f'[clang-tidy] extra content: {line}')
            self._add_issue(issue)

    def _add_issue(self, issue):
        if 'file' not in issue:
            return

        issue['file'] = os.path.relpath(issue['file'])
        if self._need_exclude(issue['file']):
            return

        formated_issue = self._format_issue(issue)

        if formated_issue['fingerprint'] in self._issues:
            return

        self._issues[formated_issue['fingerprint']] = formated_issue

    def _format_issue(self, issue_dict):
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
        issue['categories'], issue['severity'] = \
            self._get_classifiaction(type_)

        return issue


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--compiledb-path', '-p',
                    help='path to compile_commands.json', required=True)
parser.add_argument('--exclude', '-x', action='append',
                    help='path to exclude')
parser.add_argument('--clang-tidy', '-c', default='clang-tidy',
                    help='clang-tidy binary to use')
parser.add_argument('--file-list', '-f',
                    help='A file containing the list of files to check')
args, extra_args = parser.parse_known_args()

files = None  # pylint: disable=invalid-name
if args.file_list is not None:
    with open(args.file_list, 'r') as fh:
        files = fh.readlines()
        files = [file_.rstrip() for file_ in files]

formater = ClangTidy2CodeQuality(args.compiledb_path,
                                 clang_tidy=args.clang_tidy,
                                 file_list=files,
                                 excludes=args.exclude,
                                 arguments=extra_args)
formater.run()
