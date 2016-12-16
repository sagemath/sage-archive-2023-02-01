from sage.misc.superseded import deprecation
deprecation(20135, '''use "from sage.libs.gsl.all cimport *" instead of including "gsl.pxi"''')

from sage.libs.gsl.all cimport *
