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

### XXX module level docstrings in sagexc
### XXX class level docstrings in sagexc

import inspect
import os
import sagedoc

SAGE_ROOT = os.environ["SAGE_ROOT"]

def sage_getfile(obj):
    r"""
    Get the full file name associated to obj as a string.

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
    AUTHOR:
        -- Nick Alexander
    """
    source = sage_getsource(obj, is_binary=True)
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

            # Sagex often has type information; we split off the right most
            # identifier to discard this information
            s[0] = arg.split()[-1]
            # Sagex often has C pointer symbols before variable names
            s[0].lstrip('*')
            argnames.append(s[0])
            if len(s) > 1:
                argdefs.append(s[1])
        return (argnames, None, None, tuple(argdefs))
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

    AUTHOR: This is a modified version of inspect.getargspec from the
    Python Standard Library, which was taken from IPython for use in
    SAGE.
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
            return (source_lines(), 0)
        return _extract_source(source_lines, lineno)
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

def _extract_source(lines, lineno):
    r"""
    Given a list of lines or a multiline string and a starting lineno, _extract_source returns a pair ([source_lines], first_lineno).  [source_lines] is the smallest indentation block (started by 'def ', 'class ', or BOF) surrounding lineno, and first_lineno is the beginning line number of that indentation block, starting from line 0.

    EXAMPLES:
        sage: from sage.misc.sageinspect import _extract_source
        sage: s = 'import os\n# preceding comment not include\ndef test1(a, b=2):\n    if a:\n        return 1\n    return b\n# intervening comment not included\n\nclass test2():\n    pass\n    # indented comment included\n# trailing comment not included\n    \ndef test3(b,\n    a=2):\n    pass # EOF\n'
        sage: _extract_source(s, 2)
        (['def test1(a, b=2):\n', '    if a:\n', '        return 1\n', '    return b\n'], 2)
        sage: _extract_source(s, 3)
        (['def test1(a, b=2):\n', '    if a:\n', '        return 1\n', '    return b\n'], 2)
        sage: _extract_source(s, 7)
        (['def test1(a, b=2):\n', '    if a:\n', '        return 1\n', '    return b\n'], 2)
        sage: _extract_source(s, 8)
        (['class test2():\n', '    pass\n', '    # indented comment included\n'], 8)
        sage: _extract_source(s, 10)
        (['class test2():\n', '    pass\n', '    # indented comment included\n'], 8)
        sage: _extract_source(s, 13)
        (['def test3(b,\n', '    a=2):\n', '    pass # EOF\n'], 13)
    """
    if isinstance(lines, str):
        lines = lines.splitlines()
        lines = [ line + '\n' for line in lines ]

    # search backward for def or class
    # XXX won't work for module?
    # XXX breaks on embedded functions, classes, etc
    while lineno > 0 \
              and not lines[lineno].lstrip().startswith('def ') \
              and not lines[lineno].lstrip().startswith('class '): \
        lineno -= 1

    # Go down the file until we find another line indented the same
    # as the first line.
    line0 = lines[lineno]
    indent_level = len(line0) - len(line0.lstrip())
    source_lines = [line0]
    i = 1
    for line in lines[lineno+1:]:
        if len(line.strip()) == 0:
            source_lines.append(line)
            continue
        line_indent = len(line) - len(line.lstrip())
        if line_indent <= indent_level: # XXX handle strings and comments?
            break
        source_lines.append(line)

    return (source_lines, lineno)
