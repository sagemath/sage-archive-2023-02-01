"""
Support for the Notebook (introspection and setup)

AUTHORS:

- William Stein (much of this code is from IPython).
"""

import inspect
import os
import string
from cPickle import PicklingError

import sage.structure.sage_object
import sage.misc.latex
import sage.misc.pager

import sage.misc.sagedoc as sagedoc
import sage.misc.sageinspect as sageinspect

from sage.misc.preparser import preparse

import pydoc

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
    # Set this to true and plots are shown by default.
    #sage.plot.plot.SHOW_DEFAULT = True
    if object_directory:
        sage.structure.sage_object.base=object_directory
    sage.misc.latex.EMBEDDED_MODE = True
    sage.misc.pager.EMBEDDED_MODE = True

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
    Display help on s.

    .. note::

       This a wrapper around the builtin help. If formats the output
       as HTML without word wrap, which looks better in the notebook.

    INPUT:


    -  ``s`` - Python object, module, etc.


    OUTPUT: prints out help about s; it's often more more extensive
    than foo?

    TESTS::

        sage: import numpy.linalg
        sage: sage.server.support.help(numpy.linalg.norm)
        <html><table notracebacks bgcolor="#386074" cellpadding=10 cellspacing=10><tr><td bgcolor="#f5f5f5"><font color="#37546d">
        &nbsp;&nbsp;&nbsp;<a target='_new' href='cell://docs-....html'>Click to open help window</a>&nbsp;&nbsp;&nbsp;
        <br></font></tr></td></table></html>
    """
    from pydoc import resolve, html, describe
    import sage.server.notebook.interact as interact

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
    Return a list of completions in the context of globs.
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
        v = list(set(v))   # make uniq
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
    Format ``obj_name``'s docstring for printing in Sage
    notebook.

    AUTHORS:

    - William Stein: partly taken from IPython for use in Sage

    - Nick Alexander: extensioins
    """
    if system not in ['sage', 'python']:
        obj_name = system + '.' + obj_name
    try:
        obj = eval(obj_name, globs)
    except (AttributeError, NameError, SyntaxError):
        return "No object '%s' currently defined."%obj_name
    s  = ''
    try:
        filename = sageinspect.sage_getfile(obj)
        #i = filename.find('site-packages/sage/')
        #if i == -1:
        s += 'File:        %s\n'%filename
        #else:
        #    file = filename[i+len('site-packages/sage/'):]
        #    s += 'File:        <html><a href="src_browser?%s">%s</a></html>\n'%(file,file)
    except TypeError:
        pass
    s += 'Type:        %s\n'%type(obj)
    s += 'Definition:  %s\n'%sageinspect.sage_getdef(obj, obj_name)
    s += 'Docstring: \n%s\n'%sageinspect.sage_getdoc(obj, obj_name)
    return s.rstrip()

def source_code(s, globs, system='sage'):
    r"""
    Format obj's source code for printing in Sage notebook.

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
        except:
            pass
        filename = sageinspect.sage_getfile(obj)
        lines, lineno = sageinspect.sage_getsourcelines(obj, is_binary=False)
        src = ''.join(lines)
        src = sagedoc.format_src(src)
        if not lineno is None:
            src = "File: %s\nSource Code (starting at line %s):\n%s"%(filename, lineno, src)
        return src

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
    return w



def syseval(system, cmd, dir=None):
    """
    INPUT:
        system -- an object with an eval method that takes as input
                  a cmd (a string), and two dictionaries:
                           sage_globals and locals.
        dir -- an otional directory to change to before
               calling system.eval.

    OUTPUT:
        The output of system.eval is returned.

    EXAMPLES:
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
                 use_cache=False, create_local_c_file=True):
    """
    INPUT:


    -  ``filename`` - name of a file that contains cython
       code


    OUTPUT:


    -  ``module`` - the module that contains the compiled
       cython code.


    Raises an ``ImportError`` exception if anything goes
    wrong.
    """
    name, build_dir = sage.misc.cython.cython(filename, verbose=verbose,
                                            compile_message=compile_message,
                                            use_cache=use_cache,
                                            create_local_c_file=create_local_c_file)
    sys.path.append(build_dir)
    return __builtin__.__import__(name)


def cython_import_all(filename, globals, verbose=False, compile_message=False,
                     use_cache=False, create_local_c_file=True):
    """
    INPUT:


    -  ``filename`` - name of a file that contains cython
       code


    OUTPUT: changes globals using the attributes of the Cython module
    that do not begin with an underscore.

    Raises an ``ImportError`` exception if anything goes
    wrong.
    """
    m = cython_import(filename, verbose=verbose, compile_message=compile_message,
                     use_cache=use_cache,
                     create_local_c_file=create_local_c_file)
    for k, x in m.__dict__.iteritems():
        if k[0] != '_':
            globals[k] = x

