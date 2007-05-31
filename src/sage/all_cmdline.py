"""nodoctest"""

sage_mode = 'cmdline'

try:

    from sage.all import *

except ValueError, msg:
    import traceback
    t = traceback.format_exc()
    print t
    if 'type object' in str(msg):
        msg = str(msg) + '\n\n** In SAGE, the easiest fix for this problem is to type "sage -ba"\n   to rebuild all the SageX code (this takes several minutes).\n   Alternatively, touch the last .pyx file in the traceback above. **\n'
    raise ValueError, msg


from sage.calculus.predefined import *
