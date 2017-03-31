from sage.misc.superseded import deprecation
deprecation(20205, '''pari/decl.pxi is deprecated, use "from cypari2.paridecl cimport *" instead''')

from cypari2.paridecl cimport *
