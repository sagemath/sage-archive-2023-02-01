from sage.misc.superseded import deprecation
deprecation(20002, '''import "cysignals/signals.pxi" instead of "sage/ext/interrupt.pxi"''')

include "cysignals/signals.pxi"
