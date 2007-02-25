r"""
Inspect Python, Sage, and Sagex objects.

This module extends parts of python's inspect module to Sagex objects.
It does not yet handle Sagex classes or modules; a source patch to
sagexc is needed for that.

AUTHOR:
   -- originally taken from Fernando Perez's IPython
   -- modified extensively by William Stein
   -- extended by Nick Alexander
"""

import inspect
import os
import sagedoc

SAGE_ROOT = os.environ["SAGE_ROOT"]

def sage_getfile(obj):
    r"""
    Get the full file name associated to obj as a string.

    EXAMPLES:
    Works with functions, classes, and modules.
        sage: import os.path
        sage: os.path.split(sage.misc.sageinspect.sage_getfile(sage.rings.integer.Integer.factor))[1]
        'integer.pyx'
        sage: os.path.split(sage.misc.sageinspect.sage_getfile(sage.rings.integer.Integer))[1]
        'integer.pyx'
        sage: os.path.split(sage.misc.sageinspect.sage_getfile(sage.rings.integer))[1]
        'integer.pyx'
        sage: os.path.split(sage.misc.sageinspect.sage_getfile(sage.misc.sageinspect.sage_getfile))[1]
        'sageinspect.py'

    AUTHOR:
        -- Nick Alexander
    """
    try:
        # We can _extract Sagex objects' files from docstrings
        d = inspect.getdoc(obj)
        (docstring, filename, lineno) = _extract_embedded_position(d)
        return filename
    except ValueError:
        return inspect.getabsfile(obj)

def _sage_getargspec_sagex(obj):
    r"""
    inspect.getargspec for Sagex objects.

    EXAMPLES:
    If \code{obj} is not a function, we raise an exception.
        sage: sage.misc.sageinspect._sage_getargspec_sagex(sage.rings.integer.Integer)
        Traceback (most recent call last):
        ...
        TypeError: obj is not a function

    If \code{obj} is not a Sagex function, we raise an exception.
        sage: sage.misc.sageinspect._sage_getargspec_sagex(sage.misc.sageinspect.sage_getfile)
        Traceback (most recent call last):
        ...
        ValueError: Could not parse sagex argspec

    We handle default arguments, but ignore varargs and kwargs.
        sage: sage.misc.sageinspect._sage_getargspec_sagex(sage.rings.integer.Integer.factor)
        (['self', 'algorithm'], None, None, ('pari',))
        sage: sage.misc.sageinspect._sage_getargspec_sagex(sage.rings.rational.make_rational)
        (['s'], None, None, None)

    AUTHOR:
        -- Nick Alexander
    """
    if not inspect.isroutine(obj):
        raise TypeError, "obj is not a function"
    try:
        source = sage_getsource(obj, is_binary=True)
        defpos = source.find('def ')
        assert defpos > -1
        colpos = source.find(':')
        assert colpos > -1
        defsrc = source[defpos:colpos]

        lparpos = defsrc.find('(')
        assert lparpos > -1
        rparpos = defsrc.rfind(')')
        assert rparpos > -1

        argsrc = defsrc[lparpos+1:rparpos]

        # Now handle individual arguments
        # XXX this could break on embedded strings or embedded functions
        args = argsrc.split(',')

        # Now we need to take care of default arguments
        # XXX this could break on embedded strings or embedded functions with default arguments
        argnames = [] # argument names
        argdefs  = [] # default values
        for arg in args:
            s = arg.split('=')
            argname = s[0]

            # Sagex often has type information; we split off the right most
            # identifier to discard this information
            argname = argname.split()[-1]
            # Sagex often has C pointer symbols before variable names
            argname.lstrip('*')
            argnames.append(argname)
            if len(s) > 1:
                defvalue = s[1]
                # Remove quotes around strings
                defvalue = defvalue.strip('"').strip("'")
                argdefs.append(defvalue)

        if len(argdefs) > 0:
            argdefs = tuple(argdefs)
        else:
            argdefs = None

        return (argnames, None, None, argdefs)
    except:
        raise ValueError, "Could not parse sagex argspec"

def sage_getargspec(obj):
    r"""
    Return the names and default values of a function's arguments.

    A tuple of four things is returned: (args, varargs, varkw,
    defaults).  'args' is a list of the argument names (it may contain
    nested lists).  'varargs' and 'varkw' are the names of the * and
    ** arguments or None.  'defaults' is an n-tuple of the default
    values of the last n arguments.

    EXAMPLES:
        sage: sage.misc.sageinspect.sage_getargspec(sage.rings.integer.Integer.factor)
        (['self', 'algorithm'], None, None, ('pari',))
        sage: sage.misc.sageinspect.sage_getargspec(sage.misc.sageinspect.sage_getargspec)
        (['obj'], None, None, None)

    AUTHOR:
        -- William Stein: a modified version of inspect.getargspec from the
        Python Standard Library, which was taken from IPython for use in SAGE.
        -- Extensions by Nick Alexander
    """
    if inspect.isfunction(obj):
        func_obj = obj
    elif inspect.ismethod(obj):
        func_obj = obj.im_func
    else:
        # Perhaps it's binary and is defined in a Sagex file
        try:
            return _sage_getargspec_sagex(obj)
        except:
            pass
        raise TypeError, 'arg is not a Python or a Sagex function'
    args, varargs, varkw = inspect.getargs(func_obj.func_code)
    return args, varargs, varkw, func_obj.func_defaults

def sage_getdef(obj, obj_name=''):
    r"""
    Return the definition header for any callable object.

    If any exception is generated, None is returned instead and the
    exception is suppressed.

    EXAMPLES:
        sage: sage.misc.sageinspect.sage_getdef(sage.rings.rational.make_rational, obj_name='mr')
        'mr(s)'
        sage: sage.misc.sageinspect.sage_getdef(sage.rings.integer.Integer.factor, obj_name='factor')
        "factor(algorithm='pari')"

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    try:
        s = str(inspect.formatargspec(*sage_getargspec(obj)))
        s = s.strip('(').strip(')').strip()
        if s[:4] == 'self':
            s = s[4:]
        s = s.lstrip(',').strip()
        return obj_name + '(' + s + ')'
    except:
        return '%s( ... )'%obj_name

def sage_getdoc(obj, obj_name=''):
    r"""
    Return the docstring associated to obj as a string.

    If obj is a Sagex object with an embedded position in its docstring,
    the embedded position is stripped.

    EXAMPLES:
        sage: sage.misc.sageinspect.sage_getdoc(sage.rings.rational.make_rational)
        ''
        sage: sage.misc.sageinspect.sage_getdoc(sage.rings.rational).strip().splitlines()[0]
        'Rational Numbers'

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    try:
        s = sagedoc.format(str(obj._sage_doc_()))
    except AttributeError:
        s = sagedoc.format(str(obj.__doc__))
    # If there is a Sagex embedded position, it needs to be stripped
    try:
        (docstring, filename, lineno) = _extract_embedded_position(s)
        s = docstring
    except ValueError:
        pass
    if obj_name != '':
        i = obj_name.find('.')
        if i != -1:
            obj_name = obj_name[:i]
        s = s.replace('self.','%s.'%obj_name)
    return s

def sage_getsource(obj, is_binary=False):
    r"""
    Return the source code associated to obj as a string, or None.

    sage: sage.misc.sageinspect.sage_getsource(sage.rings.rational.make_rational, True)
    'def make_rational(s):\n    r = Rational()\n    r._reduce_set(s)\n    return r\n'

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    t = sage_getsourcelines(obj, is_binary)
    if not t:
        return None
    (source_lines, lineno) = t
    return ''.join(source_lines)

def sage_getsourcelines(obj, is_binary=False):
    r"""
    Return a pair ([source_lines], starting line number) of the source code associated to obj, or None.

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    if not is_binary:
        # First try the python inspection library
        try:
            return inspect.getsourcelines(obj)
        except Exception:
            pass
    # Perhaps it's binary and is defined in a Sagex file, and maybe we
    # didn't know it, so we set is_binary=False.
    try:
        try:
            d = inspect.getdoc(obj)
        except:
            return None
        (orig, filename, lineno) = _extract_embedded_position(d)
        source_lines = open(filename).readlines()
        # XXX Whole file for modules -- fails at this time because sagex does not embed positions for modules, nor for classes
        if inspect.ismodule(obj):
            return (source_lines, 0)
        else:
            return _extract_source(source_lines, lineno), lineno
    except (IOError, IndexError):
        return None

def _extract_embedded_position(docstring):
    r"""
    If docstring has a Sagex embedded position, return a tuple (original_docstring, filename, line).  If not, raise ValueError.

    EXAMPLES:
    We cannot test the filename since it depends on SAGE_ROOT.
        sage: from sage.misc.sageinspect import _extract_embedded_position
        sage: s = 'File: sage/rings/rational.pyx (starting at line 1080)\noriginal'
        sage: p = _extract_embedded_position(s)
        sage: (p[0], p[2])
        ('original', 1080)

        sage: s = 'no embedded position'
        sage: _extract_embedded_position(s)
        Traceback (most recent call last):
        ...
        ValueError: Docstring (='''no embedded position''') does not contain an embedded position

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    try:
        # isolate first line
        i = docstring.find('\n')
        if i >= 0:
            line0 = docstring[:i+1]
            original_docstring = docstring[i+1:]
        if i < 0:
            original_docstring = ''
            line0 = docstring

        if line0.startswith("File:") and '.pyx' in line0:
            # Sagex embedded position found
            parts = line0.split()
            filename = parts[1]
            assert filename.endswith('.pyx')
            last = parts[-1]
            assert last.endswith(')')
            line = last[:-1]
            filename = '%s/local/lib/python/site-packages/%s'%(SAGE_ROOT, filename)
            return (original_docstring, filename, eval(line))
    except:
        pass
    raise ValueError, "Docstring (='''%s''') does not contain an embedded position"%docstring

__test_string1 = '''
import os                                  # 1
# preceding comment not include            # 2
def test1(a, b=2):                         # 3
    if a:                                  # 4
        return 1                           # 5
    return b                               # 6
# intervening comment not included         # 7
class test2():                             # 8
    pass                                   # 9
    # indented comment not included        # 10
# trailing comment not included            # 11
def test3(b,                               # 12
          a=2):                            # 13
    pass # EOF                             # 14
'''

def _extract_source(lines, lineno):
    r"""
    Given a list of lines or a multiline string and a starting lineno,
    _extract_source returns [source_lines].  [source_lines] is the smallest
    indentation block starting at lineno.

    EXAMPLES:
        sage: from sage.misc.sageinspect import _extract_source
        sage: s = sage.misc.sageinspect.__test_string1.lstrip()
        sage: es = lambda ls, l: ''.join(_extract_source(ls, l)).rstrip()

        sage: print es(s, 3)
        def test1(a, b=2):                         # 3
            if a:                                  # 4
                return 1                           # 5
            return b                               # 6

        sage: print es(s, 8)
        class test2():                             # 8
            pass                                   # 9

        sage: print es(s, 12)
        def test3(b,                               # 12
                  a=2):                            # 13
            pass # EOF                             # 14
    """
    if lineno < 1:
        raise ValueError, "Line numbering starts at 1! (tried to extract line %s)" % lineno
    lineno -= 1

    if isinstance(lines, str):
        lines = lines.splitlines()
        lines = [ line + '\n' for line in lines ]

    return inspect.getblock(lines[lineno:])
