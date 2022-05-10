#!/usr/bin/env python3

import re
import subprocess

def run_git_command(args):
    cmd = ['git'] + args
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout = p.communicate()[0].strip().decode()
    if p.returncode != 0:
        if verbose:
            print("unable to run %s (error)" % dispcmd)
            print("stdout was %s" % stdout)
        return None, p.returncode
    return stdout, p.returncode


git_describe, rc = run_git_command(["describe", "--tags", "--dirty",
                                    "--always", "--match", "v*"])

pieces = {}
if "g" in git_describe:
    # TAG-DISTANCE-gHEX
    describe_mo = re.search(r'^(?P<tag>.+)'
                            r'-(?P<distance>\d+)'
                            r'-g(?P<short>[0-9a-f]+)'
                            r'(-(?P<dirty>dirty))?$',
                            git_describe)
    pieces['tag'] = describe_mo.group('tag')
    # distance: number of commits since tag
    pieces["distance"] = int(describe_mo.group('distance'))

    # commit: short hex revision ID
    pieces["short"] = describe_mo.group('short')
    pieces["dirty"] = describe_mo.group('dirty')
else:
    # remove prefix
    pieces['tag'] = git_describe

# major.minor.patch-prerelease+build
semver_mo = re.search(
    r'^v(?P<major>0|[1-9]\d*)'
    r'(\.(?P<minor>0|[1-9]\d*))?'
    r'(\.(?P<patch>0|[1-9]\d*))?'
    r'(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?'
    r'(?:\+(?P<build>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$',
    pieces['tag'])

if semver_mo:
    for p in ['major', 'minor', 'patch', 'prerelease', 'build']:
        if semver_mo.group(p):
            pieces[p] = semver_mo.group(p)

semver_build = []
if 'build' in pieces:
    semver_build = [pieces['build']]

if 'distance' in pieces:
    semver_build.extend([str(pieces['distance']), 'g' + pieces['short']])
    if pieces['dirty']:
        semver_build.append(pieces['dirty'])

if semver_build:
    pieces['build_part'] = '+' + '.'.join(semver_build)
else:
    pieces['build_part'] = ''

if 'prerelease' in pieces:
    pieces['prerelease'] = '-' + pieces['prerelease']
else:
    pieces['prerelease'] = ''
    
semver = '{major}.{minor}.{patch}{prerelease}{build_part}'.format(**pieces)
print(semver)
