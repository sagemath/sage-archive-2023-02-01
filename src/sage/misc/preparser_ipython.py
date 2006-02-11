###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import sage.misc.interpreter
import preparser

interface_name  = 'sage'
interface       = None
num_lines       = 0

def switch_interface(name, verbose=True):
    I = __import__('sage.interfaces.all',{},{},name)
    if not name in I.interfaces:
        raise RuntimeError, "Invalid interface %s"%name
    sage.misc.interpreter.PROMPT = '%s: '%name

    global interface, interface_name
    interface_name = name
    if name == 'sage':
        interface = None
        if verbose:
            print "\n  --> Exiting back to SAGE <-- \n"
    else:
        interface = I.__dict__[name]
        if verbose:
            print "\n  --> Switching to %s <-- \n"%interface
        if name in ['kash']:
            interface('0')
    sage.misc.interpreter.set_sage_prompt('%s: '%interface_name)

def preparse_ipython(line, reset=True):
    global num_lines

    L = line.lstrip()
    if L.startswith('%'):
        # This is an ipython magic command
        L = L[1:].strip()
        import sage.interfaces.all
        if L.lower() in sage.interfaces.all.interfaces:
            switch_interface(L.lower())
            return ''
        else:
            # only preparse non-magic lines
            return line

    if interface is None:
        return preparser.preparse(line, reset=reset)

    line = line.rstrip()
    ends_in_backslash = line.endswith('\\')
    if ends_in_backslash:
        line = line.rstrip('\\')
        num_lines += 1
    else:
        if interface_name in ['gap', 'magma', 'kash', 'singular']:
            line += ';'

    t = interface._eval_line(line,
                wait_for_prompt=(not ends_in_backslash))

    import sage.misc.interpreter
    if ends_in_backslash:

        sage.misc.interpreter.set_sage_prompt('.'*len(interface_name) + ' ')
        return ''

    else:
        sage.misc.interpreter.set_sage_prompt('%s: '%interface_name)

        # find where the real output begins
        t = '\n'.join(t.split('\n')[num_lines:])
        num_lines = 0
        #if interface_name == 'mathematica':
        #    i = t.find('Out[')
        #    t = t[i:]
        print t
        return ''
        #return 'print """%s"""'%t

