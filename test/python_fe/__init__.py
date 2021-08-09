
""" __init__.py: Python module to test the fe_engine"""

__author__ = "Nicolas Richart"
__credits__ = [
    "Nicolas Richart <nicolas.richart@epfl.ch>",
]
__copyright__ = "Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale" \
                " de Lausanne) Laboratory (LSMS - Laboratoire de Simulation" \
                " en Mécanique des Solides)"
__license__ = "LGPLv3"


import sys as __aka_sys


def export(definition):
    """
    Decorator to export definitions from sub-modules to the top-level package

    :param definition: definition to be exported
    :return: definition
    """
    __module = __aka_sys.modules[definition.__module__]
    __pkg = __aka_sys.modules[__module.__package__]
    __pkg.__dict__[definition.__name__] = definition

    if '__all__' not in __pkg.__dict__:
        __pkg.__dict__['__all__'] = []

    __pkg.__all__.append(definition.__name__)

    return definition


try:
    from termcolor import colored
except ImportError:
    # noinspection PyUnusedLocal
    def colored(string, *args, **kwargs):
        return string
__all__ = ['colored']


from . import truss_fe        # NOQA: E402
from . import static_solver   # NOQA: E402
from . import dynamic_solver  # NOQA: E402
