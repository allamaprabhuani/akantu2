#!/usr/bin/env sh

# This script requires jq, it converts the code quality json file to an output usable as errors in emacs `M-x compile`

jq -r 'unique_by(.fingerprint) | map(select(.location.lines.begin).line = .location.lines.begin) | map(select(.location.positions.begin).line = .location.positions.begin.line) | map({"path":.location.path, "line": .line, "description": .description, "engine": .engine_name, "check": .check_name}) | sort_by(.path, .line) | .[] | "/home/richart/akantu-eigen/" + (.path) + ":" + (.line | tostring) + ": " + (.description) + " [" + (.engine) + ":" + (.check) + "]"' $1
