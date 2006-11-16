"""
Support for the Notebook (introspection and setup)

AUTHOR:
   * William Stein (much of this code is from IPython).
"""

import inspect
import sage.misc.pyrex_inspect
import os
import string

import sage.plot.plot
import sage.structure.sage_object
import sage.misc.latex
import sage.misc.pager

import sage.misc.sagedoc as sagedoc


######################################################################
# Initialization
######################################################################
EMBEDDED_MODE = False
sage_globals = None
globals_at_init = None
global_names_at_init = None

def init(object_directory=None, globs={}):
    """
    Initialize SAGE for use with the web notebook interface.
    """
    global sage_globals, globals_at_init, global_names_at_init
    global EMBEDDED_MODE

    sage_globals = globs
    #globals_at_init = set(globs.keys())
    globals_at_init = globs.values()
    global_names_at_init = set(globs.keys())
    EMBEDDED_MODE = True
    sage.plot.plot.EMBEDDED_MODE = True
    # Set this to true and plots are shown by default.
    #sage.plot.plot.SHOW_DEFAULT = True
    if object_directory:
        sage.structure.sage_object.base=object_directory
    sage.misc.latex.EMBEDDED_MODE = True
    sage.misc.pager.EMBEDDED_MODE = True

    # Turn on latex print mode by default.
    #sage.misc.latex.lprint()


######################################################################
# Introspection
######################################################################

def get_rightmost_identifier(s):
    X = string.ascii_letters + string.digits + '._'
    i = len(s)-1
    while i >= 0 and s[i] in X:
        i -= 1
    return s[i+1:]

def completions(s, globs, format=False, width=90):
    """
    Return a list of completions in the context of globs.
    """
    n = len(s)
    if n == 0:
        return '(empty string)'
    if not '.' in s:
        v = [x for x in globs.keys() if x[:n] == s]
    else:
        i = s.rfind('.')
        method = s[i+1:]
        obj = s[:i]
        n = len(method)
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
        except:
            v = []
    v = list(set(v))   # make uniq
    v.sort()
    if format:
        if len(v) == 0:
            return "no completions of %s"%s
        else:
            return tabulate(v, width)
    return v

def get_argspec(obj):
    """
    Return the names and default values of a function's arguments.

    A tuple of four things is returned: (args, varargs, varkw,
    defaults).  'args' is a list of the argument names (it may contain
    nested lists).  'varargs' and 'varkw' are the names of the * and
    ** arguments or None.  'defaults' is an n-tuple of the default
    values of the last n arguments.

    AUTHOR: This is a modified version of inspect.getargspec from the
    Python Standard Library, which was taken from IPython for use in
    SAGE.
    """
    if inspect.isfunction(obj):
        func_obj = obj
    elif inspect.ismethod(obj):
        func_obj = obj.im_func
    else:
        raise TypeError, 'arg is not a Python function'
    args, varargs, varkw = inspect.getargs(func_obj.func_code)
    return args, varargs, varkw, func_obj.func_defaults


def get_def(obj, obj_name=''):
    """
    Return the definition header for any callable object.

    If any exception is generated, None is returned instead and the
    exception is suppressed.

    AUTHOR: William Stein (but taken from IPython for use in SAGE).
    """
    try:
        s = str(inspect.formatargspec(*get_argspec(obj)))
        s = s.strip('(').strip(')').strip()
        if s[:4] == 'self':
            s = s[4:]
        s = s.lstrip(',').strip()
        return obj_name + '(' + s + ')'
    except:
        return '%s( ... )'%obj_name


def get_doc(obj, obj_name=''):
    try:
        s = sagedoc.format(str(obj._sage_doc_()))
    except AttributeError:
        s = sagedoc.format(str(obj.__doc__))
    if obj_name != '':
        i = obj_name.find('.')
        if i != -1:
            obj_name = obj_name[:i]
        s = s.replace('self.','%s.'%obj_name)
    return s

def docstring(obj_name, globs):
    try:
        obj = eval(obj_name, globs)
    except (AttributeError, NameError, SyntaxError):
        return "No object '%s' currently defined."%obj_name
    s  = ''
    try:
        s += 'File:        %s\n'%inspect.getabsfile(obj)
    except TypeError:
        pass
    s += 'Type:        %s\n'%type(obj)
    s += 'Definition:  %s\n'%get_def(obj, obj_name)
    s += 'Docstring:\n%s\n'%get_doc(obj, obj_name)
    return s

def source_code(s, globs):
    """
        AUTHOR: William Stein (but partly taken from IPython for use in SAGE).
    """
    try:
        obj = eval(s, globs)
    except NameError:
        return "No object %s"%s

    try:
        try:
            fname = inspect.getabsfile(obj)
            lines, num = inspect.getsourcelines(obj)
            src = ''.join(lines)
        except TypeError:
            src = sage.misc.pyrex_inspect.getsource(obj, True)
            num =None
        src = sagedoc.format_src(src)
        if not num is None:
            src = """File: %s
Source Code (starting at line %s):\n%s"""%(fname, num, src)
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
        except (IOError, TypeError):
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



def syseval(system, cmd):
    try:
        return system.eval(cmd, locals = sage_globals)
    except TypeError:
        return system.eval(cmd)


######################################################################
# Pyrex
######################################################################
import sage.misc.pyrex
import sys
import __builtin__

def pyrex_import(filename, verbose=False, compile_message=False,
                 use_cache=False):
    """
    INPUT:
        filename -- name of a file that contains pyrex code

    OUTPUT:
        module -- the module that contains the compiled pyrex code.

    Raises an ImportError exception if anything goes wrong.
    """

    name, build_dir = sage.misc.pyrex.pyrex(filename, verbose=verbose,
                                            compile_message=compile_message,
                                            use_cache=use_cache)

    sys.path.append(build_dir)
    return __builtin__.__import__(name)


def pyrex_import_all(filename, globals, verbose=False, compile_message=False,
                     use_cache=False):
    """
    INPUT:
        filename -- name of a file that contains pyrex code

    OUTPUT:
        changes globals using the attributes of the Pyrex module
        that do not begin with an underscore.

    Raises an ImportError exception if anything goes wrong.
    """
    m = pyrex_import(filename, verbose=verbose, compile_message=compile_message,
                     use_cache=use_cache)
    for k, x in m.__dict__.iteritems():
        if k[0] != '_':
            globals[k] = x

