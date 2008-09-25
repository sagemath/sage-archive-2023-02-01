r"""
Inspect Python, Sage, and Cython objects.

This module extends parts of Python's inspect module to Cython objects.

AUTHOR:
   -- originally taken from Fernando Perez's IPython
   -- modified extensively by William Stein
   -- extended by Nick Alexander
   -- testing by Nick Alexander

EXAMPLES:
    sage: from sage.misc.sageinspect import *

Test introspection of modules defined in Python and Cython files:

    Cython modules:
        sage: sage_getfile(sage.rings.rational)
        '.../rational.pyx'

        sage: sage_getdoc(sage.rings.rational).lstrip()
        'Rational Numbers...'

        sage: sage_getsource(sage.rings.rational)[5:]
        'Rational Numbers...'

    Python modules:
        sage: sage_getfile(sage.misc.sageinspect)
        '.../sageinspect.py'

        sage: print sage_getdoc(sage.misc.sageinspect).lstrip()[:40]
        Inspect Python, Sage, and Cython objects

        sage: sage_getsource(sage.misc.sageinspect).lstrip()[5:-1]
        'Inspect Python, Sage, and Cython objects...'

Test introspection of classes defined in Python and Cython files:

    Cython classes:
        sage: sage_getfile(sage.rings.rational.Rational)
        '.../rational.pyx'

        sage: sage_getdoc(sage.rings.rational.Rational).lstrip()
        'A Rational number...'

        sage: sage_getsource(sage.rings.rational.Rational)
        'cdef class Rational...'

    Python classes:
        sage: sage_getfile(sage.misc.attach.Attach)
        '.../attach.py'

        sage: sage_getdoc(sage.misc.attach.Attach).lstrip()
        "Attach a file to a running instance of Sage..."

        sage: sage_getsource(sage.misc.attach.Attach)
        'class Attach:...'

    Python classes with no docstring, but an __init__ docstring:
        sage: class Foo:
        ...     def __init__(self):
        ...         'docstring'
        ...         pass
        ...
        sage: sage_getdoc(Foo)
        'docstring'

Test introspection of functions defined in Python and Cython files:

    Cython functions:
        sage: sage_getdef(sage.rings.rational.make_rational, obj_name='mr')
        'mr(s)'

        sage: sage_getfile(sage.rings.rational.make_rational)
        '.../rational.pyx'

        sage: sage_getdoc(sage.rings.rational.make_rational).lstrip()
        "Make a rational number ...

        sage: sage_getsource(sage.rings.rational.make_rational, True)
        'def make_rational(s):...'

    Python functions:
        sage: sage_getdef(sage.misc.sageinspect.sage_getfile, obj_name='sage_getfile')
        'sage_getfile(obj)'

        sage: sage_getfile(sage.misc.sageinspect.sage_getfile)
        '.../sageinspect.py'

        sage: sage_getdoc(sage.misc.sageinspect.sage_getfile).lstrip()
        'Get the full file name associated to obj as a string...'

        sage: sage_getsource(sage.misc.sageinspect.sage_getfile)
        'def sage_getfile(obj):...'

    Unfortunately, there is no argspec extractable from builtins:
        sage: sage_getdef(''.find, 'find')
        'find( [noargspec] )'

        sage: sage_getdef(str.find, 'find')
        'find( [noargspec] )'


"""

import inspect
import os

def isclassinstance(obj):
    r"""
    Checks if argument is instance of non built-in class
    """
    return (hasattr(obj, '__class__') and \
            hasattr(obj.__class__, '__module__') and \
            obj.__class__.__module__ not in ('__builtin__', 'exceptions'))


SAGE_ROOT = os.environ["SAGE_ROOT"]

import re
# Parse strings of form "File: sage/rings/rational.pyx (starting at line 1080)"
# "\ " protects a space in re.VERBOSE mode.
__embedded_position_re = re.compile(r'''
\A                                          # anchor to the beginning of the string
File:\ (?P<FILENAME>.*?)                    # match File: then filename
\ \(starting\ at\ line\ (?P<LINENO>\d+)\)   # match line number
\n?                                         # if there is a newline, eat it
(?P<ORIGINAL>.*)                            # the original docstring is the end
\Z                                          # anchor to the end of the string
''', re.MULTILINE | re.DOTALL | re.VERBOSE)

def _extract_embedded_position(docstring):
    r"""
    If docstring has a Cython embedded position, return a tuple (original_docstring, filename, line).  If not, return None.

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    if docstring is None:
        return None
    res = __embedded_position_re.match(docstring)
    if res is not None:
        #filename = '%s/local/lib/python/site-packages/%s' % (SAGE_ROOT, res.group('FILENAME'))
        filename = '%s/devel/sage/%s' % (SAGE_ROOT, res.group('FILENAME'))
        lineno = int(res.group('LINENO'))
        original = res.group('ORIGINAL')
        return (original, filename, lineno)
    return None

def _extract_source(lines, lineno):
    r"""
    Given a list of lines or a multiline string and a starting lineno,
    _extract_source returns [source_lines].  [source_lines] is the smallest
    indentation block starting at lineno.
    """
    if lineno < 1:
        raise ValueError, "Line numbering starts at 1! (tried to extract line %s)" % lineno
    lineno -= 1

    if isinstance(lines, str):
        lines = lines.splitlines(True) # true keeps the '\n'
    if len(lines) > 0:
        # Fixes an issue with getblock
        lines[-1] += '\n'

    return inspect.getblock(lines[lineno:])

def _sage_getargspec_cython(source):
    r"""
    inspect.getargspec from source code.

    AUTHOR:
        -- Nick Alexander
    """
    try:
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

            # Cython often has type information; we split off the right most
            # identifier to discard this information
            argname = argname.split()[-1]
            # Cython often has C pointer symbols before variable names
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
        raise ValueError, "Could not parse cython argspec"

def sage_getfile(obj):
    r"""
    Get the full file name associated to obj as a string.

    AUTHOR:
        -- Nick Alexander
    """
    # We try to extract from docstrings, because Python's inspect
    # will happily report compiled .so files
    d = inspect.getdoc(obj)
    pos = _extract_embedded_position(d)
    if pos is not None:
        (_, filename, _) = pos
        return filename

    # The instance case
    if isclassinstance(obj):
        return inspect.getabsfile(obj.__class__)
    # No go? fall back to inspect.
    return inspect.getabsfile(obj)

def sage_getargspec(obj):
    r"""
    Return the names and default values of a function's arguments.

    A tuple of four things is returned: (args, varargs, varkw,
    defaults).  'args' is a list of the argument names (it may contain
    nested lists).  'varargs' and 'varkw' are the names of the * and
    ** arguments or None.  'defaults' is an n-tuple of the default
    values of the last n arguments.

    AUTHOR:
        -- William Stein: a modified version of inspect.getargspec from the
        Python Standard Library, which was taken from IPython for use in SAGE.
        -- Extensions by Nick Alexander
    """
    if not callable(obj):
        raise TypeError, "obj is not a code object"

    if inspect.isfunction(obj):
        func_obj = obj
    elif inspect.ismethod(obj):
        func_obj = obj.im_func
    elif isclassinstance(obj):
        return sage_getargspec(obj.__class__.__call__)
    elif inspect.isclass(obj):
        return sage_getargspec(obj.__call__)
    else:
        # Perhaps it is binary and defined in a Cython file
        source = sage_getsource(obj, is_binary=True)
        if source:
            return _sage_getargspec_cython(source)
        else:
            func_obj = obj

    # Otherwise we're (hopefully!) plain Python, so use inspect
    try:
        args, varargs, varkw = inspect.getargs(func_obj.func_code)
    except AttributeError:
        args, varargs, varkw = inspect.getargs(func_obj)
    try:
        defaults = func_obj.func_defaults
    except AttributeError:
        defaults = tuple([])
    return args, varargs, varkw,

def sage_getdef(obj, obj_name=''):
    r"""
    Return the definition header for any callable object.

    If an exception is generated, None is returned instead and the
    exception is suppressed.

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    try:
        spec = sage_getargspec(obj)
        s = str(inspect.formatargspec(*spec))
        s = s.strip('(').strip(')').strip()
        if s[:4] == 'self':
            s = s[4:]
        s = s.lstrip(',').strip()
        return obj_name + '(' + s + ')'
    except (AttributeError, TypeError, ValueError):
        return '%s( [noargspec] )'%obj_name

def sage_getdoc(obj, obj_name=''):
    r"""
    Return the docstring associated to obj as a string.

    If obj is a Cython object with an embedded position in its docstring,
    the embedded position is stripped.

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    if obj is None: return ''
    import sage.misc.sagedoc
    r = None
    try:
        r = obj._sage_doc_()
    except AttributeError:
        r = obj.__doc__

    #Check to see if there is an __init__ method, and if there
    #is, use its docstring.
    if r is None and hasattr(obj, '__init__'):
        r = obj.__init__.__doc__

    if r is None:
        return ''
    s = sage.misc.sagedoc.format(str(r))

    # If there is a Cython embedded position, it needs to be stripped
    pos = _extract_embedded_position(s)
    if pos is not None:
        s, _, _ = pos

    # Fix object naming
    if obj_name != '':
        i = obj_name.find('.')
        if i != -1:
            obj_name = obj_name[:i]
        s = s.replace('self.','%s.'%obj_name)

    return s

def sage_getsource(obj, is_binary=False):
    r"""
    Return the source code associated to obj as a string, or None.

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
    Return a pair ([source_lines], starting line number) of the source
    code associated to obj, or None.

    At this time we ignore is_binary in favour of a 'do our best' strategy.

    AUTHOR:
        -- William Stein
        -- Extensions by Nick Alexander
    """
    # Check if we deal with instance
    if isclassinstance(obj):
        obj=obj.__class__
    # If we can handle it, we do.  This is because Python's inspect will
    # happily dump binary for cython extension source code.
    d = inspect.getdoc(obj)
    pos = _extract_embedded_position(d)
    if pos is None:
        return inspect.getsourcelines(obj)

    (orig, filename, lineno) = pos
    try:
        source_lines = open(filename).readlines()
    except IOError:
        return None

    return _extract_source(source_lines, lineno), lineno



__internal_teststring = '''
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
    pass # EOF                             # 14'''

def __internal_tests():
    r"""
    Test internals of the sageinspect module.

    sage: from sage.misc.sageinspect import *
    sage: from sage.misc.sageinspect import _extract_source, _extract_embedded_position, _sage_getargspec_cython, __internal_teststring

    If docstring is None, nothing bad happens:
        sage: sage_getdoc(None)
        ''

        sage: sage_getsource(sage)
        "...all..."

    A cython function with default arguments:
        sage: sage_getdef(sage.rings.integer.Integer.factor, obj_name='factor')
        "factor(algorithm='pari', proof='True', limit='None')"

    A cython method without an embedded position can lead to surprising errors:
        sage: sage_getsource(sage.rings.integer.Integer.__init__, is_binary=True)
        Traceback (most recent call last):
        ...
        TypeError: arg is not a module, class, method, function, traceback, frame, or code object

        sage: sage_getdef(sage.rings.integer.Integer.__init__, obj_name='__init__')
        '__init__( [noargspec] )'

    Test _extract_source with some likely configurations, including no trailing
    newline at the end of the file:

        sage: s = __internal_teststring.strip()
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

    Test _sage_getargspec_cython with multiple default arguments and a type:
        sage: _sage_getargspec_cython("def init(self, x=None, base=0):")
        (['self', 'x', 'base'], None, None, ('None', '0'))
        sage: _sage_getargspec_cython("def __init__(self, x=None, base=0):")
        (['self', 'x', 'base'], None, None, ('None', '0'))
        sage: _sage_getargspec_cython("def __init__(self, x=None, unsigned int base=0):")
        (['self', 'x', 'base'], None, None, ('None', '0'))

    Test _extract_embedded_position:
        We cannot test the filename since it depends on SAGE_ROOT.

        Make sure things work with no trailing newline:
            >>> _extract_embedded_position('File: sage/rings/rational.pyx (starting at line 1080)')
            ('', '.../rational.pyx', 1080)

        And with a trailing newline:
            >>> s = 'File: sage/rings/rational.pyx (starting at line 1080)\n'
            >>> _extract_embedded_position(s)
            ('', '.../rational.pyx', 1080)

        And with an original docstring:
            >>> s = 'File: sage/rings/rational.pyx (starting at line 1080)\noriginal'
            >>> _extract_embedded_position(s)
            ('original', '.../rational.pyx', 1080)

        And with a complicated original docstring:
            >>> s = 'File: sage/rings/rational.pyx (starting at line 1080)\n\n\noriginal test\noriginal'
            >>> _extract_embedded_position(s)
            ('\n\noriginal test\noriginal', ..., 1080)

            >>> s = 'no embedded position'
            >>> _extract_embedded_position(s) is None
            True
    """
