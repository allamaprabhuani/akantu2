__copyright__ = "Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"


__all__ = []

from .py_engine import *  # noqa: F401,F403
from . import py_engine as __pye  # noqa: F401

__all__.append(__pye.__all__)
