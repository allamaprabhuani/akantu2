#!/usr/bin/env python3

""" test_patch_linear_heat_transfer_explicit.py: heat transfer patch test in
python"""

__author__ = "Guillaume Anciaux"
__credits__ = [
    "Guillaume Anciaux <guillaume.anciaux@epfl.ch>",
]
__copyright__ = "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"

import sys
from patch_test_linear_heat_transfer_fixture import TestPatchTestHTMLinear
from patch_test_linear_heat_transfer_fixture import run_test_generic
import akantu


def run_test(self_):
    run_test_generic(self_, akantu._explicit_lumped_mass)


def test():
    TestPatchTestHTMLinear.TYPED_TEST(run_test, "Explicit")


if 'pytest' not in sys.modules:
    test()
