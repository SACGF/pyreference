from __future__ import absolute_import

from .gene import *
from .mirna import *
from .reference import *
from .referenceargparse import *
from .transcript import *

__version__ = "0.6.3"


def get_json_schema_version():
    """ Return an int which increments upon breaking changes - ie anything other than patch """
    major, minor, patch = __version__.split(".")
    return 1000 * int(major) + int(minor)
