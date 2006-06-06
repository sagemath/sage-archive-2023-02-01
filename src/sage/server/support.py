import inspect
import os


import sage.plot.plot
import sage.ext.sage_object
import sage.misc.latex

import sage.misc.sagedoc as sagedoc


######################################################################
# Initialization
######################################################################

def init(object_directory):
    """
    Initialize SAGE for use with the web notebook interface.
    """
    sage.plot.plot.EMBEDDED_MODE = True
    sage.ext.sage_object.base=object_directory
    sage.misc.latex.EMBEDDED_MODE = True



######################################################################
# Introspection
######################################################################

def completions(s, globs, format=False, width=90):
    """
    Return a list of completions in the context of globs.
    """
    n = len(s)
    if not '.' in s:
        v = [x for x in globs.keys() if x[:n] == s]
    else:
        i = s.rfind('.')
        method = s[i+1:]
        obj = s[:i]
        n = len(method)
        if globs.has_key(obj):
            O = globs[obj]
            D = dir(O)
            try:
                D += O.trait_names()
            except (AttributeError, TypeError):
                pass
            if method == '':
                v = [obj + '.'+x for x in D if x and x[0] != '_']
            else:
                v = [obj + '.'+x for x in D if x[:n] == method]
        else:
            v = []
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

    AUTHOR: Taken from IPython for use in SAGE.
    """
    try:
        return obj_name + inspect.formatargspec(*get_argspec(obj))
    except:
        return None


def get_doc(obj):
    try:
        return sagedoc.format(str(obj._sage_doc_()))
    except AttributeError:
        return sagedoc.format(str(obj.__doc__))

def docstring(obj_name, globs):
    try:
        obj = eval(obj_name, globs)
    except (AttributeError, NameError, SyntaxError):
        return "No object '%s' currently defined."%obj_name
    s  = 'Type:        %s\n'%type(obj)
    s += 'Definition:  %s\n'%get_def(obj, obj_name)
    s += 'Docstring:\n%s\n'%get_doc(obj)
    return s

def source_code(s, globs):
    obj = eval(s, globs)
    try:
        z = inspect.getsourcelines(obj)[0]
    except TypeError:
        return "Source code not available."
    return '\n'.join(z)


def tabulate(v, width=90):
    if len(v) == 0:
        return ''
    n = max(len(x) for x in v)
    ncols = width // (n + 7)
    nrows = len(v)//ncols + 1
    col_width = width//ncols
    s = ''
    for r in range(nrows):
        for c in range(ncols):
            i = r + c*nrows
            if i < len(v):
                w = v[i]
                s += w + ' '*(col_width - len(w))
        s += '\n'
    return s
