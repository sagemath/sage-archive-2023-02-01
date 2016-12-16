from sage.misc.superseded import deprecation
deprecation(20158, '''python.pxi is deprecated, use "from cpython cimport *" instead''')

from cpython.ref cimport *
from cpython.exc cimport *
from cpython.module cimport *
from cpython.mem cimport *
from cpython.tuple cimport *
from cpython.list cimport *
from cpython.object cimport *
from cpython.sequence cimport *
from cpython.mapping cimport *
from cpython.iterator cimport *
from cpython.number cimport *
from cpython.int cimport *
from cpython.bool cimport *
from cpython.long cimport *
from cpython.float cimport *
from cpython.complex cimport *
from cpython.string cimport *
from cpython.dict cimport *
from cpython.instance cimport *
from cpython.function cimport *
from cpython.method cimport *
from cpython.set cimport *
from cpython.slice cimport *
