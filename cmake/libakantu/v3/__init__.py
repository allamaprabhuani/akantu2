
""" __init__.py: module for gdb pretty printer for akantu v3"""

__author__ = "Nicolas Richart"
__credits__ = [
    "Nicolas Richart <nicolas.richart@epfl.ch>",
]
__copyright__ = "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"


import gdb

# Load the pretty-printers.
from .printers import register_akantu_printers
register_akantu_printers(gdb.current_objfile())
