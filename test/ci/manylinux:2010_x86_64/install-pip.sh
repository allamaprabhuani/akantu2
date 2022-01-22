#!/usr/bin/env bash
set -euo pipefail

for PIP in /opt/python/*/bin/pip; do
  $PIP install --disable-pip-version-check --upgrade pip
  $PIP install scikit-build==0.8.1
done
