#!/usr/bin/env python3
import sys
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

from .issue_generator_clang_tidy import ClangTidyIssueGenerator
from .issue_generator_clang_format import ClangFormatIssueGenerator
from .issue_generator_warnings import WarningsIssueGenerator
from .issue_list import IssueList

def run(cmd, **kwargs):
    import json
    issue_list = IssueList(**kwargs)

    if cmd == 'clang_tidy':
        tool = ClangTidyIssueGenerator(issue_list, **kwargs)
    elif cmd == 'clang_format':
        tool = ClangFormatIssueGenerator(issue_list, **kwargs)
    elif cmd == 'warnings':
        tool = WarningsIssueGenerator(issue_list, **kwargs)

    tool.generate_issues()

    print(json.dumps(issue_list.issues))
