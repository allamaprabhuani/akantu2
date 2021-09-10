#!/usr/bin/env python3

from . import print_debug, print_info
from .issue_generator import IssueGenerator
import os
import re
import copy
import json
import subprocess


class ClangToolIssueGenerator(IssueGenerator):
    """issue generator for clang tidy"""

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

    def __init__(self, issues_list, tool, **kwargs):
        self._tool = tool
        opts = copy.copy(kwargs)
        super().__init__(issues_list, **kwargs)

        compiledb_path = opts.pop('compiledb_path')
        arguments = opts.pop('arguments', None)
        clang_tool = opts.pop('clang_tool_executable', tool)
        file_list = opts.pop('file_list', None)

        self._command = [clang_tool]

        if 'need_compiledb' in opts and opts['need_compiledb']:
            self._command.extend(['-p', compiledb_path])
        if arguments is not None:
            self._command.extend(arguments)

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
            need_exclude = self._issues._need_exclude(filename)
            if need_exclude:
                print_debug(f'[{self._tool}] exluding file: {filename}')
                continue
            print_info(f'[{self._tool}] adding file: {filename}')
            self._files.append(filename)

    def _run_command(self, command):
        print_info(f'''[{self._tool}] command: {' '.join(command)}''')
        popen = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.DEVNULL,
                                 universal_newlines=True)

        for stdout_line in iter(popen.stdout.readline, ""):
            clean_line = self.ANSI_ESCAPE.sub('', stdout_line)
            yield clean_line

        popen.stdout.close()

        return_code = popen.wait()
        if return_code:
            print_debug(
                f"[{self._tool}] {command} ReturnCode {return_code}")
