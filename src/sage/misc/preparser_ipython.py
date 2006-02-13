###########################################################################
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

import sage.misc.interpreter

import log
import preparser

interface_name  = 'sage'
interface       = None
num_lines       = 0

# TODO: Do this right?
magma_colon_equals = False

q_lines = []

def switch_interface(name, verbose=True, control_d=False):
    I = __import__('sage.interfaces.all',{},{},name)
    if not name in I.interfaces:
        raise RuntimeError, "Invalid interface %s"%name
    sage.misc.interpreter.PROMPT = '%s'%name

    global interface, interface_name
    interface_name = name
    if name == 'sage':
        interface = None
        if control_d:
            log.offset -= 1
        if verbose:
            print "\n  --> Exiting back to SAGE <-- \n"
    else:
        if not (interface is None):
            interface._post_interact()
        interface = I.__dict__[name]
        if verbose:
            print "\n  --> Switching to %s <-- \n"%interface
        interface._pre_interact()
        if name in ['kash']:
            interface('0')
    sage.misc.interpreter.set_sage_prompt('%s'%interface_name)

def preparse_ipython(line, reset=True):
    global num_lines
    global q_lines

    L = line.lstrip()
    if L.startswith('%'):
        # This should be installed as an Ipython magic command,
        # but I don't know how yet...
        L = L[1:].strip()
        import sage.interfaces.all
        if L.lower() in sage.interfaces.all.interfaces:
            switch_interface(L.lower())
            return "''"
        else:
            # only preparse non-magic lines
            return line

    if interface is None:
        return preparser.preparse(line, reset=reset)

    line = preparse_imports_from_sage(interface, line)
    line = line.rstrip()
    ends_in_backslash = line.endswith('\\')
    if ends_in_backslash:
        line = line.rstrip('\\')
        num_lines += 1
    else:
        if interface_name in ['gap', 'magma', 'kash', 'singular']:
            line += ';'
            if magma_colon_equals and interface_name == 'magma':
                line = line.replace(':=','=').replace('=',':=')

    if ends_in_backslash:
        q_lines.append(line)
    else:
        if len(q_lines) > 0:
            line = ''.join(q_lines) + line
        q_lines = []
        # TODO: do sage substitutions here
        #t = interface._eval_line(line)
        t = interface.eval(line)

    import sage.misc.interpreter
    if ends_in_backslash:

        sage.misc.interpreter.set_sage_prompt('.'*len(interface_name))

    else:

        sage.misc.interpreter.set_sage_prompt('%s'%interface_name)
        #print t

    return 'r"""%s"""'%t


_v_ = None

def preparse_imports_from_sage(interface, line, locals={}):
    """
    The input line is being fed to the given interface.
    This function extracts any "sage(zzz)"'s, parses
    them, and replaces them by appropriate objects
    in the interface.   This is used for moving
    objects from SAGE back into the interface.
    """
    import sage_eval
    i = line.find('%s('%interface_name)
    if i == -1:
        return line
    j = i + line[i:].find(')')
    expr = line[i+len(interface_name)+1:j]
    expr = preparser.preparse(expr)
    s = 'import sage.misc.preparser_ipython; \
         sage.misc.preparser_ipython._v_ = %s(%s)'%(interface_name, expr)
    #print s
    __IPYTHON__.runsource(s)
    #print _v_
    line = line[:i] + _v_.name() + line[j+2:]
    return line

