"""
Support for Notebook Introspection and Setup

AUTHORS:

- William Stein (much of this code is from IPython).

- Nick Alexander
"""

import os
import string
from cPickle import PicklingError

import sage.structure.sage_object
import sage.misc.latex
import sage.misc.pager

import sage.misc.sagedoc as sagedoc
import sage.misc.sageinspect as sageinspect

from sage.misc.preparser import preparse

######################################################################
# Initialization
######################################################################
EMBEDDED_MODE = False
sage_globals = None
globals_at_init = None
global_names_at_init = None

def init(object_directory=None, globs={}):
    r"""
    Initialize Sage for use with the web notebook interface.
    """
    global sage_globals, globals_at_init, global_names_at_init
    global EMBEDDED_MODE

    os.environ['PAGER'] = 'cat'

    sage_globals = globs
    #globals_at_init = set(globs.keys())
    globals_at_init = globs.values()
    global_names_at_init = set(globs.keys())
    EMBEDDED_MODE = True

    import sage.plot.plot
    sage.plot.plot.EMBEDDED_MODE = True
    if object_directory:
        sage.structure.sage_object.base=object_directory
    sage.misc.latex.EMBEDDED_MODE = True
    sage.misc.pager.EMBEDDED_MODE = True
    sage.misc.sageinspect.EMBEDDED_MODE = True

    setup_systems(globs)
    sage.misc.session.init(globs)


def setup_systems(globs):
    from sage.misc.inline_fortran import InlineFortran
    fortran = InlineFortran(globs)
    globs['fortran'] = fortran


######################################################################
# Introspection
######################################################################
def help(obj):
    """
    Display HTML help for ``obj``, a Python object, module, etc.  This
    help is often more extensive than that given by 'obj?'.  This
    function does not return a value --- it prints HTML as a side
    effect.

    .. note::

       This a wrapper around the built-in help. If formats the output
       as HTML without word wrap, which looks better in the notebook.

    INPUT:

    -  ``obj`` - a Python object, module, etc.

    TESTS::

        sage: import numpy.linalg
        sage: import os, sage.misc.misc ; current_dir = os.getcwd()
        sage: os.chdir(sage.misc.misc.tmp_dir('server_doctest'))
        sage: sage.server.support.help(numpy.linalg.norm)
        <html><table notracebacks bgcolor="#386074" cellpadding=10 cellspacing=10><tr><td bgcolor="#f5f5f5"><font color="#37546d">
        &nbsp;&nbsp;&nbsp;<a target='_new' href='cell://docs-....html'>Click to open help window</a>&nbsp;&nbsp;&nbsp;
        <br></font></tr></td></table></html>
        sage: os.chdir(current_dir)
    """
    from pydoc import resolve, html, describe

    print '<html><table notracebacks bgcolor="#386074" cellpadding=10 cellspacing=10><tr><td bgcolor="#f5f5f5"><font color="#37546d">'
    object, name = resolve(obj)
    page = html.page(describe(object), html.document(object, name))
    page = page.replace('<a href','<a ')
    n = 0
    while True:
        filename = 'docs-%s.html'%n
        if not os.path.exists(filename): break
        n += 1
    open(filename, 'w').write(page)
    print "&nbsp;&nbsp;&nbsp;<a target='_new' href='cell://%s'>Click to open help window</a>&nbsp;&nbsp;&nbsp;"%filename
    print '<br></font></tr></td></table></html>'

def get_rightmost_identifier(s):
    X = string.ascii_letters + string.digits + '._'
    i = len(s)-1
    while i >= 0 and s[i] in X:
        i -= 1
    return s[i+1:]

def completions(s, globs, format=False, width=90, system="None"):
    """
    Return a list of completions in the given context.

    INPUT:

    - ``globs`` - a string:object dictionary; context in which to
      search for completions, e.g., :func:`globals()`

    - ``format`` - a bool (default: False); whether to tabulate the
      list

    - ``width`` - an int; character width of the table

    - ``system`` - a string (default: 'None'); system prefix for the
      completions
    """
    if system not in ['sage', 'python']:
        prepend = system + '.'
        s = prepend + s
    else:
        prepend = ''
    n = len(s)
    if n == 0:
        return '(empty string)'
    try:
        if not '.' in s and not '(' in s:
            v = [x for x in globs.keys() if x[:n] == s] + \
                [x for x in __builtins__.keys() if x[:n] == s]
        else:
            if not ')' in s:
                i = s.rfind('.')
                method = s[i+1:]
                obj = s[:i]
                n = len(method)
            else:
                obj = preparse(s)
                method = ''
            try:
                O = eval(obj, globs)
                D = dir(O)
                try:
                    D += O.trait_names()
                except (AttributeError, TypeError):
                    pass
                if method == '':
                    v = [obj + '.'+x for x in D if x and x[0] != '_']
                else:
                    v = [obj + '.'+x for x in D if x[:n] == method]
            except Exception, msg:
                v = []
        v = list(set(v))   # make unique
        v.sort()
    except Exception, msg:
        v = []

    if prepend:
        i = len(prepend)
        v = [x[i:] for x in v]

    if format:
        if len(v) == 0:
            return "No completions of '%s' currently defined"%s
        else:
            return tabulate(v, width)
    return v

def docstring(obj_name, globs, system='sage'):
    r"""
    Format an object's docstring to process and display in the Sage
    notebook.

    INPUT:

    - ``obj_name`` - a string; a name of an object

    - ``globs`` - a string:object dictionary; a context in which to
      evaluate ``obj_name``

    - ``system`` - a string (default: 'sage'); the system to which to
      confine the search

    OUTPUT:

    - a string containing the object's file, type, definition, and
      docstring or a message stating the object is not defined

    AUTHORS:

    - William Stein: partly taken from IPython for use in Sage

    - Nick Alexander: extensions
    """
    if system not in ['sage', 'python']:
        obj_name = system + '.' + obj_name
    try:
        obj = eval(obj_name, globs)
    except (AttributeError, NameError, SyntaxError):
        return "No object '%s' currently defined."%obj_name
    s  = ''
    newline = "\n\n"  # blank line to start new paragraph
    try:
        filename = sageinspect.sage_getfile(obj)
        #i = filename.find('site-packages/sage/')
        #if i == -1:
        s += '**File:** %s'%filename
        s += newline
        #else:
        #    file = filename[i+len('site-packages/sage/'):]
        #    s += 'File:        <html><a href="src_browser?%s">%s</a></html>\n'%(file,file)
    except TypeError:
        pass
    s += '**Type:** %s'%type(obj)
    s += newline
    s += '**Definition:** %s'%sageinspect.sage_getdef(obj, obj_name)
    s += newline
    s += '**Docstring:**'
    s += newline
    s += sageinspect.sage_getdoc(obj, obj_name)
    return s.rstrip()

def source_code(s, globs, system='sage'):
    r"""
    Format an object's source code to process and display in the
    Sage notebook.

    INPUT:

    - ``s`` - a string; a name of an object

    - ``globs`` - a string:object dictionary; a context in which to
      evaluate ``s``

    - ``system`` - a string (default: 'sage'); the system to which to
      confine the search

    OUTPUT:

    - a string containing the object's file, starting line number, and
      source code

    AUTHORS:

    - William Stein: partly taken from IPython for use in Sage

    - Nick Alexander: extensions
    """
    if system not in ['sage', 'python']:
        s = system + '.' + s

    try:
        obj = eval(s, globs)
    except NameError:
        return "No object %s"%s

    try:
        try:
            return obj._sage_src_()
        except Exception:
            pass
        newline = "\n\n"  # blank line to start new paragraph
        indent = "    "   # indent source code to mark it as a code block

        filename = sageinspect.sage_getfile(obj)
        lines, lineno = sageinspect.sage_getsourcelines(obj, is_binary=False)
        src = indent.join(lines)
        src = indent + sagedoc.format_src(src)
        if not lineno is None:
            output = "**File:** %s"%filename
            output += newline
            output += "**Source Code** (starting at line %s)::"%lineno
            output += newline
            output += src
        return output

    except (TypeError, IndexError), msg:
        print msg
        return "Source code for %s not available."%obj

def tabulate(v, width=90, ncols=3):
    e = len(v)
    if e == 0:
        return ''
    while True:
        col_widths = []
        nrows = e//ncols + 1
        for c in range(ncols):
            m = max([0] + [len(v[r+c*nrows]) for r in range(nrows) if r+c*nrows < e])
            col_widths.append(m+3)
        if ncols > 1 and max(col_widths + [0]) > width//ncols:
            ncols -= 1
        else:
            break
    n = max(len(x) for x in v)
    s = ''
    for r in range(nrows):
        for c in range(ncols):
            i = r + c*nrows
            if i < e:
                w = v[i]
                s += w + ' '*(col_widths[c] - len(w))
        s += '\n'
    return s

def save_session(filename):
    D = {}
    v = variables(with_types=False)
    for k in v:
        x = sage_globals[k]
        try:
            _ = sage.structure.sage_object.loads(sage.structure.sage_object.dumps(x))
        except (IOError, TypeError, PicklingError):
            if k != 'fortran':  # this is a hack to get around the inline fortran object being
                                # *incredibly* hackish in how it is implemented; the right
                                # fix is to rewrite the fortran inline to *not* be so incredibly
                                # hackish.  See trac #2891.
                print "Unable to save %s"%k
        else:
            D[k] = x
    print "Saving variables to object %s.sobj"%filename
    sage.structure.sage_object.save(D, filename)

def load_session(v, filename, state):
    D = {}
    for k, x in v.iteritems():
        try:
            _ = sage.structure.sage_object.loads(sage.structure.sage_object.dumps(x))
        except (IOError, TypeError):
            print "Unable to save %s"%k
        else:
            D[k] = x
    print "Saving variables to %s"%filename
    sage.structure.sage_object.save(D, filename)

def _is_new_var(x, v):
    if x[:2] == '__':
        return False
    if not x in global_names_at_init:
        return True

    # You might think this would take a long time
    # since globals_at_init has several thousand entries.
    # However, it takes 0.0 seconds, which is not noticeable
    # given that there is at least 0.1 seconds delay
    # when refreshing the web page!
    for y in globals_at_init:
        if v is y:
            return False
    return True

def variables(with_types=True):
    if with_types:
        w = ['%s-%s'%(x,type(v)) for x, v in sage_globals.iteritems() if \
             _is_new_var(x, v)]
    else:
        w = [x for x, v in sage_globals.iteritems() if \
             _is_new_var(x, v)]
    w.sort()
    return tuple(w)



def syseval(system, cmd, dir=None):
    """
    Evaluate an input with a "system" object that can evaluate inputs
    (e.g., python, gap).

    INPUT:

    - ``system`` - an object with an eval method that takes an input

    - ``cmd`` - a string input

    - ``sage_globals`` - a string:object dictionary

    - dir - a string (default: None); an optional directory to change
      to before calling :func:`system.eval`

    OUTPUT:

    - :func:`system.eval`'s output

    EXAMPLES::

        sage: from sage.misc.python import python
        sage: sage.server.support.syseval(python, '2+4/3')
        3
        ''
        sage: sage.server.support.syseval(python, 'import os; os.chdir(".")')
        ''
        sage: sage.server.support.syseval(python, 'import os; os.chdir(1,2,3)')
        Traceback (most recent call last):
        ...
        TypeError: chdir() takes exactly 1 argument (3 given)
        sage: sage.server.support.syseval(gap, "2+3")
        '5'
    """
    if dir:
        if hasattr(system.__class__, 'chdir'):
            system.chdir(dir)
    return system.eval(cmd, sage_globals, locals = sage_globals)

######################################################################
# Cython
######################################################################
import sage.misc.cython
import sys
import __builtin__

def cython_import(filename, verbose=False, compile_message=False,
                 use_cache=False, create_local_c_file=True, **kwds):
    """
    Compile a file containing Cython code, then import and return the
    module.  Raises an ``ImportError`` if anything goes wrong.

    INPUT:

    - ``filename`` - a string; name of a file that contains Cython
      code

    See the function :func:`sage.misc.cython.cython` for documentation
    for the other inputs.

    OUTPUT:

    - the module that contains the compiled Cython code.
    """
    name, build_dir = sage.misc.cython.cython(filename, verbose=verbose,
                                            compile_message=compile_message,
                                            use_cache=use_cache,
                                            create_local_c_file=create_local_c_file,
                                            **kwds)
    sys.path.append(build_dir)
    return __builtin__.__import__(name)


def cython_import_all(filename, globals, verbose=False, compile_message=False,
                     use_cache=False, create_local_c_file=True):
    """
    Imports all non-private (i.e., not beginning with an underscore)
    attributes of the specified Cython module into the given context.
    This is similar to::

        from module import *

    Raises an ``ImportError`` exception if anything goes wrong.

    INPUT:

    - ``filename`` - a string; name of a file that contains Cython
      code
    """
    m = cython_import(filename, verbose=verbose, compile_message=compile_message,
                     use_cache=use_cache,
                     create_local_c_file=create_local_c_file)
    for k, x in m.__dict__.iteritems():
        if k[0] != '_':
            globals[k] = x

