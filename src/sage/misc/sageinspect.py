r"""
Inspect Python, Sage, and Cython objects.

This module extends parts of Python's inspect module to Cython objects.

AUTHORS:

- originally taken from Fernando Perez's IPython
- William Stein (extensive modifications)
- Nick Alexander (extensions)
- Nick Alexander (testing)

EXAMPLES::

    sage: from sage.misc.sageinspect import *

Test introspection of modules defined in Python and Cython files:

Cython modules::

    sage: sage_getfile(sage.rings.rational)
    '.../rational.pyx'

    sage: sage_getdoc(sage.rings.rational).lstrip()
    'Rational Numbers...'

    sage: sage_getsource(sage.rings.rational)[5:]
    'Rational Numbers...'

Python modules::

    sage: sage_getfile(sage.misc.sageinspect)
    '.../sageinspect.py'

    sage: print sage_getdoc(sage.misc.sageinspect).lstrip()[:40]
    Inspect Python, Sage, and Cython objects

    sage: sage_getsource(sage.misc.sageinspect).lstrip()[5:-1]
    'Inspect Python, Sage, and Cython objects...'

Test introspection of classes defined in Python and Cython files:

Cython classes::

    sage: sage_getfile(sage.rings.rational.Rational)
    '.../rational.pyx'

    sage: sage_getdoc(sage.rings.rational.Rational).lstrip()
    'A Rational number...'

    sage: sage_getsource(sage.rings.rational.Rational)
    'cdef class Rational...'

Python classes::

    sage: import sage.misc.attach
    sage: sage_getfile(sage.misc.attach.Attach)
    '.../attach.py'

    sage: sage_getdoc(sage.misc.attach.Attach).lstrip()
    'Attach a file to a running instance of Sage...'

    sage: sage_getsource(sage.misc.attach.Attach)
    'class Attach:...'

Python classes with no docstring, but an __init__ docstring::

    sage: class Foo:
    ...     def __init__(self):
    ...         'docstring'
    ...         pass
    ...
    sage: sage_getdoc(Foo)
    'docstring\n'

Test introspection of functions defined in Python and Cython files:

Cython functions::

    sage: sage_getdef(sage.rings.rational.make_rational, obj_name='mr')
    'mr(s)'

    sage: sage_getfile(sage.rings.rational.make_rational)
    '.../rational.pyx'

    sage: sage_getdoc(sage.rings.rational.make_rational).lstrip()
    "Make a rational number ..."

    sage: sage_getsource(sage.rings.rational.make_rational, True)[4:]
    'make_rational(s):...'

Python functions::

    sage: sage_getdef(sage.misc.sageinspect.sage_getfile, obj_name='sage_getfile')
    'sage_getfile(obj)'

    sage: sage_getfile(sage.misc.sageinspect.sage_getfile)
    '.../sageinspect.py'

    sage: sage_getdoc(sage.misc.sageinspect.sage_getfile).lstrip()
    "Get the full file name associated to ``obj`` as a string..."

    sage: sage_getsource(sage.misc.sageinspect.sage_getfile)[4:]
    'sage_getfile(obj):...'

Unfortunately, there is no argspec extractable from builtins::

    sage: sage_getdef(''.find, 'find')
    'find( [noargspec] )'

    sage: sage_getdef(str.find, 'find')
    'find( [noargspec] )'
"""

import ast
import inspect
import os
import tokenize
EMBEDDED_MODE = False

def isclassinstance(obj):
    r"""
    Checks if argument is instance of non built-in class

    INPUT: ``obj`` - object

    EXAMPLES::

        sage: from sage.misc.sageinspect import isclassinstance
        sage: isclassinstance(int)
        False
        sage: isclassinstance(FreeModule)
        True
        sage: isclassinstance(SteenrodAlgebra)
        True
        sage: class myclass: pass
        sage: isclassinstance(myclass)
        False
        sage: class mymetaclass(type): pass
        sage: class myclass2:
        ...       __metaclass__ = mymetaclass
        sage: isclassinstance(myclass2)
        False
    """
    return (not inspect.isclass(obj) and \
            hasattr(obj, '__class__') and \
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
    If docstring has a Cython embedded position, return a tuple
    (original_docstring, filename, line).  If not, return None.

    INPUT: ``docstring`` (string)

    EXAMPLES::

       sage: from sage.misc.sageinspect import _extract_embedded_position
       sage: import inspect
       sage: _extract_embedded_position(inspect.getdoc(var))[1][-21:]
       'sage/calculus/var.pyx'

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


class BlockFinder:
    """
    Provide a tokeneater() method to detect the end of a code block.

    This is the Python library's :class:`inspect.BlockFinder` modified
    to recognize Cython definitions.
    """
    def __init__(self):
        self.indent = 0
        self.islambda = False
        self.started = False
        self.passline = False
        self.last = 1

    def tokeneater(self, type, token, srow_scol, erow_ecol, line):
        srow, scol = srow_scol
        erow, ecol = erow_ecol
        if not self.started:
            # look for the first "(cp)def", "class" or "lambda"
            if token in ("def", "cpdef", "class", "lambda"):
                if token == "lambda":
                    self.islambda = True
                self.started = True
            self.passline = True    # skip to the end of the line
        elif type == tokenize.NEWLINE:
            self.passline = False   # stop skipping when a NEWLINE is seen
            self.last = srow
            if self.islambda:       # lambdas always end at the first NEWLINE
                raise inspect.EndOfBlock
        elif self.passline:
            pass
        elif type == tokenize.INDENT:
            self.indent = self.indent + 1
            self.passline = True
        elif type == tokenize.DEDENT:
            self.indent = self.indent - 1
            # the end of matching indent/dedent pairs end a block
            # (note that this only works for "def"/"class" blocks,
            #  not e.g. for "if: else:" or "try: finally:" blocks)
            if self.indent <= 0:
                raise inspect.EndOfBlock
        elif self.indent == 0 and type not in (tokenize.COMMENT, tokenize.NL):
            # any other token on the same indentation level end the previous
            # block as well, except the pseudo-tokens COMMENT and NL.
            raise inspect.EndOfBlock

def _getblock(lines):
    """
    Extract the block of code at the top of the given list of lines.

    This is the Python library's :func:`inspect.getblock`, except that
    it uses an instance of our custom :class:`BlockFinder`.
    """
    blockfinder = BlockFinder()
    try:
        tokenize.tokenize(iter(lines).next, blockfinder.tokeneater)
    except (inspect.EndOfBlock, IndentationError):
        pass
    return lines[:blockfinder.last]

def _extract_source(lines, lineno):
    r"""
    Given a list of lines or a multiline string and a starting lineno,
    _extract_source returns [source_lines].  [source_lines] is the smallest
    indentation block starting at lineno.

    INPUT:

    - ``lines`` - string or list of strings
    - ``lineno`` - positive integer

    EXAMPLES::

        sage: from sage.misc.sageinspect import _extract_source
        sage: s2 = "#hello\n\n  class f():\n    pass\n\n#goodbye"
        sage: _extract_source(s2, 3)
        ['  class f():\n', '    pass\n']
    """
    if lineno < 1:
        raise ValueError, "Line numbering starts at 1! (tried to extract line %s)" % lineno
    lineno -= 1

    if isinstance(lines, str):
        lines = lines.splitlines(True) # true keeps the '\n'
    if len(lines) > 0:
        # Fixes an issue with getblock
        lines[-1] += '\n'

    return _getblock(lines[lineno:])


class SageArgSpecVisitor(ast.NodeVisitor):
    """
    A simple visitor class that walks an abstract-syntax tree (AST)
    for a Python function's argspec.  It returns the contents of nodes
    representing the basic Python types: None, booleans, numbers,
    strings, lists, tuples, and dictionaries.  We use this class in
    :func:`_sage_getargspec_from_ast` to extract an argspec from a
    function's or method's source code.

    EXAMPLES::

        sage: import ast, sage.misc.sageinspect as sms
        sage: visitor = sms.SageArgSpecVisitor()
        sage: visitor.visit(ast.parse('[1,2,3]').body[0].value)
        [1, 2, 3]
        sage: visitor.visit(ast.parse("{'a':('e',2,[None,({False:True},'pi')]), 37.0:'temp'}").body[0].value)
        {'a': ('e', 2, [None, ({False: True}, 'pi')]), 37.0: 'temp'}
        sage: v = ast.parse("jc = ['veni', 'vidi', 'vici']").body[0]; v
        <_ast.Assign object at ...>
        sage: [x for x in dir(v) if not x.startswith('__')]
        ['_attributes', '_fields', 'col_offset', 'lineno', 'targets', 'value']
        sage: visitor.visit(v.targets[0])
        'jc'
        sage: visitor.visit(v.value)
        ['veni', 'vidi', 'vici']
    """
    def visit_Name(self, node):
        """
        Visit a Python AST :class:`ast.Name` node.

        INPUT:

        - ``node`` - the node instance to visit

        OUTPUT:

        - None, True, False, or the ``node``'s name as a string.

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Name(ast.parse(x).body[0].value)
            sage: [vis(n) for n in ['True', 'False', 'None', 'foo', 'bar']]
            [True, False, None, 'foo', 'bar']
            sage: [type(vis(n)) for n in ['True', 'False', 'None', 'foo', 'bar']]
            [<type 'bool'>, <type 'bool'>, <type 'NoneType'>, <type 'str'>, <type 'str'>]
        """
        what = node.id
        if what == 'None':
            return None
        elif what == 'True':
            return True
        elif what == 'False':
            return False
        return node.id

    def visit_Num(self, node):
        """
        Visit a Python AST :class:`ast.Num` node.

        INPUT:

        - ``node`` - the node instance to visit

        OUTPUT:

        - the number the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Num(ast.parse(x).body[0].value)
            sage: [vis(n) for n in ['123', '0.0', str(-pi.n())]]
            [123, 0.0, -3.14159265358979]
        """
        return node.n

    def visit_Str(self, node):
        r"""
        Visit a Python AST :class:`ast.Str` node.

        INPUT:

        - ``node`` - the node instance to visit

        OUTPUT:

        - the string the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Str(ast.parse(x).body[0].value)
            sage: [vis(s) for s in ['"abstract"', "u'syntax'", '''r"tr\ee"''']]
            ['abstract', u'syntax', 'tr\\ee']
        """
        return node.s

    def visit_List(self, node):
        """
        Visit a Python AST :class:`ast.List` node.

        INPUT:

        - ``node`` - the node instance to visit

        OUTPUT:

        - the list the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_List(ast.parse(x).body[0].value)
            sage: [vis(l) for l in ['[]', "['s', 't', 'u']", '[[e], [], [pi]]']]
            [[], ['s', 't', 'u'], [['e'], [], ['pi']]]
         """
        t = []
        for n in node.elts:
            t.append(self.visit(n))
        return t

    def visit_Tuple(self, node):
        """
        Visit a Python AST :class:`ast.Tuple` node.

        INPUT:

        - ``node`` - the node instance to visit

        OUTPUT:

        - the tuple the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Tuple(ast.parse(x).body[0].value)
            sage: [vis(t) for t in ['()', '(x,y)', '("Au", "Al", "Cu")']]
            [(), ('x', 'y'), ('Au', 'Al', 'Cu')]
        """
        t = []
        for n in node.elts:
            t.append(self.visit(n))
        return tuple(t)

    def visit_Dict(self, node):
        """
        Visit a Python AST :class:`ast.Dict` node.

        INPUT:

        - ``node`` - the node instance to visit

        OUTPUT:

        - the dictionary the ``node`` represents

        EXAMPLES::

            sage: import ast, sage.misc.sageinspect as sms
            sage: visitor = sms.SageArgSpecVisitor()
            sage: vis = lambda x: visitor.visit_Dict(ast.parse(x).body[0].value)
            sage: [vis(d) for d in ['{}', "{1:one, 'two':2, other:bother}"]]
            [{}, {1: 'one', 'other': 'bother', 'two': 2}]
        """
        d = {}
        for k, v in zip(node.keys, node.values):
            d[self.visit(k)] = self.visit(v)
        return d

def _sage_getargspec_from_ast(source):
    r"""
    Return an argspec for a Python function or method by compiling its
    source to an abstract-syntax tree (AST) and walking its ``args``
    subtrees with :class:`SageArgSpecVisitor`.  We use this in
    :func:`_sage_getargspec_cython`.

    INPUT:

    - ``source`` - a string; the function's (or method's) source code
      definition.  The function's body is ignored.

    OUTPUT:

    - an instance of :obj:`inspect.ArgSpec`, i.e., a named tuple

    EXAMPLES::

        sage: import inspect, sage.misc.sageinspect as sms
        sage: from_ast = sms._sage_getargspec_from_ast
        sage: s = "def f(a, b=2, c={'a': [4, 5.5, False]}, d=(None, True)):\n    return"
        sage: from_ast(s)
        ArgSpec(args=['a', 'b', 'c', 'd'], varargs=None, keywords=None, defaults=(2, {'a': [4, 5.5, False]}, (None, True)))
        sage: context = {}
        sage: exec compile(s, '<string>', 'single') in context
        sage: inspect.getargspec(context['f'])
        ArgSpec(args=['a', 'b', 'c', 'd'], varargs=None, keywords=None, defaults=(2, {'a': [4, 5.5, False]}, (None, True)))
        sage: from_ast(s) == inspect.getargspec(context['f'])
        True
        sage: set(from_ast(sms.sage_getsource(x)) == inspect.getargspec(x) for x in [factor, identity_matrix, Graph.__init__])
        set([True])
    """
    ast_args = ast.parse(source.lstrip()).body[0].args

    visitor = SageArgSpecVisitor()
    args = [visitor.visit(a) for a in ast_args.args]
    defaults = [visitor.visit(d) for d in ast_args.defaults]

    return inspect.ArgSpec(args, ast_args.vararg, ast_args.kwarg,
                           tuple(defaults) if defaults else None)

def _sage_getargspec_cython(source):
    r"""
    inspect.getargspec from source code.  That is, get the names and
    default values of a function's arguments.

    INPUT: ``source`` - a string of Cython code

    OUTPUT: a tuple (``args``, None, None, ``argdefs``), where
    ``args`` is the list of arguments and ``argdefs`` is their default
    values (as strings, so 2 is represented as '2', etc.).

    EXAMPLES::

        sage: from sage.misc.sageinspect import _sage_getargspec_cython as sgc
        sage: sgc("cpdef double abc(self, x=None, base=0):")
        (['self', 'x', 'base'], None, None, (None, 0))
        sage: sgc("def __init__(self, x=None, unsigned int base=0):")
        (['self', 'x', 'base'], None, None, (None, 0))
        sage: sgc('def o(p, *q, r={}, **s) except? -1:')
        (['p', '*q', 'r', '**s'], None, None, ({},))
        sage: sgc('cpdef how(r=(None, "u:doing?")):')
        ArgSpec(args=['r'], varargs=None, keywords=None, defaults=((None, 'u:doing?'),))
        sage: sgc('def _(x="):"):')
        ArgSpec(args=['x'], varargs=None, keywords=None, defaults=('):',))
        sage: sgc('def f(z = {(1,2,3): True}):\n    return z')
        ArgSpec(args=['z'], varargs=None, keywords=None, defaults=({(1, 2, 3): True},))
        sage: sgc('def f(double x, z = {(1,2,3): True}):\n    return z')
        Traceback (most recent call last):
        ...
        ValueError: Could not parse cython argspec

    AUTHOR:

    - Nick Alexander
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
            # only process arg if it has positive length
            if len(arg) > 0:
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
                    # eval defvalue so we aren't just returning strings
                    try:
                        argdefs.append(eval(defvalue))
                    except NameError:
                        argdefs.append(defvalue)

        if len(argdefs) > 0:
            argdefs = tuple(argdefs)
        else:
            argdefs = None

        return (argnames, None, None, argdefs)

    except Exception:
        try:
            # Try to parse the entire definition as Python and get an
            # argspec.
            beg = re.search(r'def([ ]+\w+)+[ ]*\(', source).end() - 1
            proxy = 'def dummy' + source[beg:] + '\n    return'
            return _sage_getargspec_from_ast(proxy)

        except Exception:
            try:
                # Try to parse just the arguments as a Python argspec.
                beg = re.search(r'def([ ]+\w+)+[ ]*\(', source).end() - 1
                end = re.search(r'\)[ ]*:', source).end()
                proxy = 'def dummy' + source[beg:end] + '\n    return'
                return _sage_getargspec_from_ast(proxy)

            except Exception:
                raise ValueError, "Could not parse cython argspec"

def sage_getfile(obj):
    r"""
    Get the full file name associated to ``obj`` as a string.

    INPUT: ``obj``, a Sage object, module, etc.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getfile
        sage: sage_getfile(sage.rings.rational)[-23:]
        'sage/rings/rational.pyx'
        sage: sage_getfile(Sq)[-41:]
        'sage/algebras/steenrod_algebra_element.py'

    AUTHOR:
    - Nick Alexander
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

    INPUT: ``obj``, a function

    OUTPUT: A tuple of four things is returned: ``(args, varargs, varkw,
    defaults)``.  ``args`` is a list of the argument names (it may contain
    nested lists).  ``varargs`` and ``varkw`` are the names of the * and
    ** arguments or None.  ``defaults`` is an n-tuple of the default
    values of the last n arguments.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getargspec
        sage: sage_getargspec(identity_matrix)
        (['ring', 'n', 'sparse'], None, None, (0, False))
        sage: sage_getargspec(Poset)
        (['data', 'element_labels', 'cover_relations'], None, None, (None, None, False))
        sage: sage_getargspec(factor)
        (['n', 'proof', 'int_', 'algorithm', 'verbose'],
         None,
         'kwds',
         (None, False, 'pari', 0))

    AUTHORS:

    - William Stein: a modified version of inspect.getargspec from the
      Python Standard Library, which was taken from IPython for use in Sage.
    - Extensions by Nick Alexander
    """
    from sage.misc.lazy_attribute import lazy_attribute
    from sage.misc.abstract_method import AbstractMethod
    if isinstance(obj, (lazy_attribute, AbstractMethod)):
        source = sage_getsource(obj)
        return _sage_getargspec_cython(source)
    if not callable(obj):
        raise TypeError, "obj is not a code object"
    if hasattr(obj, '_sage_argspec_'):
        return obj._sage_argspec_()
    if inspect.isfunction(obj):
        func_obj = obj
    elif inspect.ismethod(obj):
        func_obj = obj.im_func
    elif isclassinstance(obj):
        return sage_getargspec(obj.__class__.__call__)
    elif inspect.isclass(obj):
        return sage_getargspec(obj.__call__)
    elif (hasattr(obj, '__objclass__') and hasattr(obj, '__name__') and
          obj.__name__ == 'next'):
        # Handle sage.rings.ring.FiniteFieldIterator.next and similar
        # slot wrappers.  This is mainly to suppress Sphinx warnings.
        return ['self'], None, None, None
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
    return args, varargs, varkw, defaults

def sage_getdef(obj, obj_name=''):
    r"""
    Return the definition header for any callable object.

    INPUT:

    - ``obj`` - function
    - ``obj_name`` - string (optional, default '')

    ``obj_name`` is prepended to the output.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getdef
        sage: sage_getdef(identity_matrix)
        '(ring, n=0, sparse=False)'
        sage: sage_getdef(identity_matrix, 'identity_matrix')
        'identity_matrix(ring, n=0, sparse=False)'

    Check that trac ticket #6848 has been fixed::

        sage: sage_getdef(RDF.random_element)
        '(min=-1, max=1)'

    If an exception is generated, None is returned instead and the
    exception is suppressed.

    AUTHORS:

    - William Stein
    - extensions by Nick Alexander
    """
    try:
        spec = sage_getargspec(obj)
        s = str(inspect.formatargspec(*spec))
        s = s.strip('(').strip(')').strip()
        if s[:4] == 'self':
            s = s[4:]
        s = s.lstrip(',').strip()
        # for use with typesetting the definition with the notebook:
        # sometimes s contains "*args" or "**keywds", and the
        # asterisks confuse ReST/sphinx/docutils, so escape them:
        # change * to \*, and change ** to \**.
        if EMBEDDED_MODE:
            s = s.replace('**', '\\**')  # replace ** with \**
            t = ''
            while True:  # replace * with \*
                i = s.find('*')
                if i == -1:
                    break
                elif i > 0 and s[i-1] == '\\':
                    if s[i+1] == "*":
                        t += s[:i+2]
                        s = s[i+2:]
                    else:
                        t += s[:i+1]
                        s = s[i+1:]
                    continue
                elif i > 0 and s[i-1] == '*':
                    t += s[:i+1]
                    s = s[i+1:]
                    continue
                else:
                    t += s[:i] + '\\*'
                    s = s[i+1:]
            s = t + s
        return obj_name + '(' + s + ')'
    except (AttributeError, TypeError, ValueError):
        return '%s( [noargspec] )'%obj_name

def _sage_getdoc_unformatted(obj):
    r"""
    Return the unformatted docstring associated to ``obj`` as a
    string.  Feed the results from this into the
    sage.misc.sagedoc.format for printing to the screen.

    INPUT: ``obj``, a function, module, etc.: something with a docstring.

    If ``obj`` is a Cython object with an embedded position in its
    docstring, the embedded position is stripped.

    EXAMPLES::

        sage: from sage.misc.sageinspect import _sage_getdoc_unformatted
        sage: _sage_getdoc_unformatted(identity_matrix)[5:44]
        'Return the `n \\times n` identity matrix'

    TESTS:

    Test that we suppress useless built-in output (Ticket #3342)

        sage: from sage.misc.sageinspect import _sage_getdoc_unformatted
        sage: _sage_getdoc_unformatted(isinstance.__class__)
        ''

    AUTHORS:

    - William Stein
    - extensions by Nick Alexander
    """
    if obj is None: return ''
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

    # Check if the __doc__ attribute was actually a string, and
    # not a 'getset_descriptor' or similar.
    import types
    if not isinstance(r, types.StringTypes):
        return ''

    from sagenb.misc.misc import encoded_str
    return encoded_str(r)

def sage_getdoc(obj, obj_name='', embedded_override=False):
    r"""
    Return the docstring associated to ``obj`` as a string.

    INPUT: ``obj``, a function, module, etc.: something with a docstring.

    If ``obj`` is a Cython object with an embedded position in its
    docstring, the embedded position is stripped.

    If optional argument ``embedded_override`` is False (its default
    value), then the string is formatted according to the value of
    EMBEDDED_MODE.  If this argument is True, then it is formatted as
    if EMBEDDED_MODE were True.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getdoc
        sage: sage_getdoc(identity_matrix)[3:39]
        'Return the n x n identity matrix ove'

    AUTHORS:

    - William Stein
    - extensions by Nick Alexander
    """
    import sage.misc.sagedoc
    if obj is None: return ''
    r = _sage_getdoc_unformatted(obj)

    if r is None:
        return ''

    s = sage.misc.sagedoc.format(str(r), embedded=(embedded_override or EMBEDDED_MODE))

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

    INPUT:

    - ``obj`` - function, etc.
    - ``is_binary`` - boolean, ignored

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getsource
        sage: sage_getsource(identity_matrix, True)[4:45]
        'identity_matrix(ring, n=0, sparse=False):'
        sage: sage_getsource(identity_matrix, False)[4:45]
        'identity_matrix(ring, n=0, sparse=False):'

    AUTHORS:

    - William Stein
    - extensions by Nick Alexander
    """
    #First we should check if the object has a _sage_src_
    #method.  If it does, we just return the output from
    #that.  This is useful for getting pexpect interface
    #elements to behave similar to regular Python objects
    #with respect to introspection.
    try:
        return obj._sage_src_()
    except (AttributeError, TypeError):
        pass

    t = sage_getsourcelines(obj, is_binary)
    if not t:
        return None
    (source_lines, lineno) = t
    return ''.join(source_lines)

def sage_getsourcelines(obj, is_binary=False):
    r"""
    Return a pair ([source_lines], starting line number) of the source
    code associated to obj, or None.

    INPUT:

    - ``obj`` - function, etc.
    - ``is_binary`` - boolean, ignored

    OUTPUT: (source_lines, lineno) or None: ``source_lines`` is a list
    of strings, and ``lineno`` is an integer.

    At this time we ignore ``is_binary`` in favour of a 'do our best' strategy.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getsourcelines
        sage: sage_getsourcelines(matrix, True)[1]
        34
        sage: sage_getsourcelines(matrix, False)[0][0][4:]
        'matrix(*args, **kwds):\n'

    AUTHORS:

    - William Stein
    - Extensions by Nick Alexander
    """

    try:
        return obj._sage_src_lines_()
    except (AttributeError, TypeError):
        pass

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

def sage_getvariablename(obj, omit_underscore_names=True):
    """
    Attempt to get the name of a Sage object.

    INPUT:

    - ``obj`` - an object
    - ``omit_underscore_names`` (optional, default True)

    If the user has assigned an object ``obj`` to a variable name,
    then return that variable name.  If several variables point to
    ``obj``, return a sorted list of those names.  If
    ``omit_underscore_names`` is True (the default) then omit names
    starting with an underscore "_".

    This is a modified version of code taken from
    http://pythonic.pocoo.org/2009/5/30/finding-objects-names,
    written by Georg Brandl.

    EXAMPLES::

        sage: from sage.misc.sageinspect import sage_getvariablename
        sage: A = random_matrix(ZZ, 100)
        sage: sage_getvariablename(A)
        'A'
        sage: B = A
        sage: sage_getvariablename(A)
        ['A', 'B']

    If an object is not assigned to a variable, an empty list is returned::

        sage: sage_getvariablename(random_matrix(ZZ, 60))
        []
    """
    import gc
    result = []
    for referrer in gc.get_referrers(obj):
        if isinstance(referrer, dict):
            for k, v in referrer.iteritems():
                if v is obj:
                    if isinstance(k, str):
                        if (not omit_underscore_names) or not k.startswith('_'):
                            result.append(k)
    if len(result) == 1:
        return result[0]
    else:
        return sorted(result)

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

    EXAMPLES::

        sage: from sage.misc.sageinspect import *
        sage: from sage.misc.sageinspect import _extract_source, _extract_embedded_position, _sage_getargspec_cython, __internal_teststring

    If docstring is None, nothing bad happens::

        sage: sage_getdoc(None)
        ''

        sage: sage_getsource(sage)
        "...all..."

    A cython function with default arguments (one of which is a string)::

        sage: sage_getdef(sage.rings.integer.Integer.factor, obj_name='factor')
        "factor(algorithm='pari', proof=True, limit=None)"

    A cython method without an embedded position can lead to surprising errors::

        sage: sage_getsource(sage.rings.integer.Integer.__init__, is_binary=True)
        Traceback (most recent call last):
        ...
        TypeError: arg is not a module, class, method, function, traceback, frame, or code object

        sage: sage_getdef(sage.rings.integer.Integer.__init__, obj_name='__init__')
        '__init__( [noargspec] )'

    Test _extract_source with some likely configurations, including no trailing
    newline at the end of the file::

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

    Test _sage_getargspec_cython with multiple default arguments and a type::

        sage: _sage_getargspec_cython("def init(self, x=None, base=0):")
        (['self', 'x', 'base'], None, None, (None, 0))
        sage: _sage_getargspec_cython("def __init__(self, x=None, base=0):")
        (['self', 'x', 'base'], None, None, (None, 0))
        sage: _sage_getargspec_cython("def __init__(self, x=None, unsigned int base=0):")
        (['self', 'x', 'base'], None, None, (None, 0))

    Test _extract_embedded_position:

    We cannot test the filename since it depends on SAGE_ROOT.

    Make sure things work with no trailing newline::

        sage: _extract_embedded_position('File: sage/rings/rational.pyx (starting at line 1080)')
        ('', '.../rational.pyx', 1080)

    And with a trailing newline::

        sage: s = 'File: sage/rings/rational.pyx (starting at line 1080)\n'
        sage: _extract_embedded_position(s)
        ('', '.../rational.pyx', 1080)

    And with an original docstring::

        sage: s = 'File: sage/rings/rational.pyx (starting at line 1080)\noriginal'
        sage: _extract_embedded_position(s)
        ('original', '.../rational.pyx', 1080)

    And with a complicated original docstring::

        sage: s = 'File: sage/rings/rational.pyx (starting at line 1080)\n\n\noriginal test\noriginal'
        sage: _extract_embedded_position(s)
        ('\n\noriginal test\noriginal', ..., 1080)

        sage: s = 'no embedded position'
        sage: _extract_embedded_position(s) is None
        True
    """
