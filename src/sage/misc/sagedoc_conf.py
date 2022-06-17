r"""
Sphinx configuration shared by sage.misc.sphinxify and sage_docbuild
"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys
import os

from docutils import nodes
from docutils.transforms import Transform
from sphinx.ext.doctest import blankline_re

# The reST default role (used for this markup: `text`) to use for all documents.
default_role = 'math'

#####################################################

def process_docstring_aliases(app, what, name, obj, options, docstringlines):
    """
    Change the docstrings for aliases to point to the original object.
    """
    basename = name.rpartition('.')[2]
    if hasattr(obj, '__name__') and obj.__name__ != basename:
        docstringlines[:] = ['See :obj:`%s`.' % name]

def process_directives(app, what, name, obj, options, docstringlines):
    """
    Remove 'nodetex' and other directives from the first line of any
    docstring where they appear.
    """
    if len(docstringlines) == 0:
        return
    first_line = docstringlines[0]
    directives = [ d.lower() for d in first_line.split(',') ]
    if 'nodetex' in directives:
        docstringlines.pop(0)

def process_docstring_cython(app, what, name, obj, options, docstringlines):
    """
    Remove Cython's filename and location embedding.
    """
    if len(docstringlines) <= 1:
        return

    first_line = docstringlines[0]
    if first_line.startswith('File:') and '(starting at' in first_line:
        #Remove the first two lines
        docstringlines.pop(0)
        docstringlines.pop(0)

def process_docstring_module_title(app, what, name, obj, options, docstringlines):
    """
    Removes the first line from the beginning of the module's docstring.  This
    corresponds to the title of the module's documentation page.
    """
    if what != "module":
        return

    #Remove any additional blank lines at the beginning
    title_removed = False
    while len(docstringlines) > 1 and not title_removed:
        if docstringlines[0].strip() != "":
            title_removed = True
        docstringlines.pop(0)

    #Remove any additional blank lines at the beginning
    while len(docstringlines) > 1:
        if docstringlines[0].strip() == "":
            docstringlines.pop(0)
        else:
            break

skip_picklability_check_modules = [
    #'sage.misc.test_nested_class', # for test only
    'sage.misc.latex',
    'sage.misc.explain_pickle',
    '__builtin__',
]

def check_nested_class_picklability(app, what, name, obj, skip, options):
    """
    Print a warning if pickling is broken for nested classes.
    """
    if hasattr(obj, '__dict__') and hasattr(obj, '__module__'):
        # Check picklability of nested classes.  Adapted from
        # sage.misc.nested_class.modify_for_nested_pickle.
        module = sys.modules[obj.__module__]
        for (nm, v) in obj.__dict__.items():
            if (isinstance(v, type) and
                v.__name__ == nm and
                v.__module__ == module.__name__ and
                getattr(module, nm, None) is not v and
                v.__module__ not in skip_picklability_check_modules):
                # OK, probably this is an *unpicklable* nested class.
                app.warn('Pickling of nested class %r is probably broken. '
                         'Please set the metaclass of the parent class to '
                         'sage.misc.nested_class.NestedClassMetaclass.' % (
                        v.__module__ + '.' + name + '.' + nm))

def skip_member(app, what, name, obj, skip, options):
    """
    To suppress Sphinx warnings / errors, we

    - Don't include [aliases of] builtins.

    - Don't include the docstring for any nested class which has been
      inserted into its module by
      :class:`sage.misc.NestedClassMetaclass` only for pickling.  The
      class will be properly documented inside its surrounding class.

    - Optionally, check whether pickling is broken for nested classes.

    - Optionally, include objects whose name begins with an underscore
      ('_'), i.e., "private" or "hidden" attributes, methods, etc.

    Otherwise, we abide by Sphinx's decision.  Note: The object
    ``obj`` is excluded (included) if this handler returns True
    (False).
    """
    if 'SAGE_CHECK_NESTED' in os.environ:
        check_nested_class_picklability(app, what, name, obj, skip, options)

    if getattr(obj, '__module__', None) == '__builtin__':
        return True

    objname = getattr(obj, "__name__", None)
    if objname is not None:
        # check if name was inserted to the module by NestedClassMetaclass
        if name.find('.') != -1 and objname.find('.') != -1:
            if objname.split('.')[-1] == name.split('.')[-1]:
                return True

    if 'SAGE_DOC_UNDERSCORE' in os.environ:
        if name.split('.')[-1].startswith('_'):
            return False

    return skip

def process_dollars(app, what, name, obj, options, docstringlines):
    r"""
    Replace dollar signs with backticks.

    See sage.misc.sagedoc.process_dollars for more information.
    """
    if len(docstringlines) and name.find("process_dollars") == -1:
        from sage.misc.sagedoc import process_dollars as sagedoc_dollars
        s = sagedoc_dollars("\n".join(docstringlines))
        lines = s.split("\n")
        for i in range(len(lines)):
            docstringlines[i] = lines[i]

def process_inherited(app, what, name, obj, options, docstringlines):
    """
    If we're including inherited members, omit their docstrings.
    """
    if not options.get('inherited-members'):
        return

    if what in ['class', 'data', 'exception', 'function', 'module']:
        return

    name = name.split('.')[-1]

    if what == 'method' and hasattr(obj, 'im_class'):
        if name in obj.im_class.__dict__.keys():
            return

    if what == 'attribute' and hasattr(obj, '__objclass__'):
        if name in obj.__objclass__.__dict__.keys():
            return

    for i in range(len(docstringlines)):
        docstringlines.pop()

def skip_TESTS_block(app, what, name, obj, options, docstringlines):
    """
    Skip blocks labeled "TESTS:".

    See sage.misc.sagedoc.skip_TESTS_block for more information.
    """
    from sage.misc.sagedoc import skip_TESTS_block as sagedoc_skip_TESTS
    if not docstringlines:
        # No docstring, so don't do anything. See Trac #19932.
        return
    s = sagedoc_skip_TESTS("\n".join(docstringlines))
    lines = s.split("\n")
    for i in range(len(lines)):
        docstringlines[i] = lines[i]
    while len(docstringlines) > len(lines):
        del docstringlines[len(lines)]

class SagemathTransform(Transform):
    """
    Transform for code-blocks.

    This allows Sphinx to treat code-blocks with prompt "sage:" as
    associated with the pycon lexer, and in particular, to change
    "<BLANKLINE>" to a blank line.
    """
    default_priority = 500

    def apply(self):
        for node in self.document.traverse(nodes.literal_block):
            if node.get('language') is None and node.astext().startswith('sage:'):
                node['language'] = 'ipycon'
                source = node.rawsource
                source = blankline_re.sub('', source)
                node.rawsource = source
                node[:] = [nodes.Text(source)]

from sage.misc.sageinspect import sage_getargspec
autodoc_builtin_argspec = sage_getargspec

# This is only used by sage.misc.sphinxify
def setup(app):
    app.connect('autodoc-process-docstring', process_docstring_cython)
    app.connect('autodoc-process-docstring', process_directives)
    app.connect('autodoc-process-docstring', process_docstring_module_title)
    app.connect('autodoc-process-docstring', process_dollars)
    app.connect('autodoc-process-docstring', process_inherited)
    if os.environ.get('SAGE_SKIP_TESTS_BLOCKS', False):
        app.connect('autodoc-process-docstring', skip_TESTS_block)
    app.connect('autodoc-skip-member', skip_member)
    app.add_transform(SagemathTransform)
