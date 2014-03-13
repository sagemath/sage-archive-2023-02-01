r"""
Some tools for developers

AUTHORS:

- Nicolas M. Thiery: initial version

- Vincent Delecroix (2012 and 2013): improve import_statements
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def module_names_cmp(x,y):
    r"""
    A comparison function for module names.

    This function first compares the depth of the modules and then
    breaks ties by alphabetical order.

    .. SEEALSO:: This function is used in :func:`import_statements`.

    TESTS::

        sage: from sage.misc.dev_tools import module_names_cmp
        sage: l = ['sage.groups.perm_gps', 'sage.combinat', 'sage.all', 'sage.plot.plot3d']
        sage: sorted(l, cmp=module_names_cmp)
        ['sage.all', 'sage.combinat', 'sage.groups.perm_gps', 'sage.plot.plot3d']
    """
    test = cmp(x.count('.'), y.count('.'))
    if test: return test
    return cmp(x,y)

def runsnake(command):
    """
    Graphical profiling with ``runsnake``

    INPUT:

    - ``command`` -- the command to be run as a string.

    EXAMPLES::

        sage: runsnake("list(SymmetricGroup(3))")        # optional - runsnake

    ``command`` is first preparsed (see :func:`preparse`)::

        sage: runsnake('for x in range(1,4): print x^2') # optional - runsnake
        1
        4
        9

    :func:`runsnake` requires the program ``runsnake``. Due to non
    trivial dependencies (python-wxgtk, ...), installing it within the
    Sage distribution is unpractical. Hence, we recommend installing
    it with the system wide Python. On Ubuntu 10.10, this can be done
    with::

        > sudo apt-get install python-profiler python-wxgtk2.8 python-setuptools
        > sudo easy_install RunSnakeRun

    See the ``runsnake`` website for instructions for other platforms.

    :func:`runsnake` further assumes that the system wide Python is
    installed in ``/usr/bin/python``.

    .. seealso::

        - `The runsnake website <http://www.vrplumber.com/programming/runsnakerun/>`_
        - ``%prun``
        - :class:`Profiler`

    """
    import cProfile, os
    from sage.misc.misc import tmp_filename, get_main_globals
    from sage.misc.preparser import preparse
    tmpfile = tmp_filename()
    cProfile.runctx(preparse(command.lstrip().rstrip()), get_main_globals(), locals(), filename=tmpfile)
    os.system("/usr/bin/python -E `which runsnake` %s &"%tmpfile)

def print_or_update(string, data):
    r"""
    if ``data`` is ``None`` then print the string otherwise append ``string`` to
    the data.
    """
    if data is None:
        print string
    else:
        data.append(string)

def print_import_statement(module, name, lazy, answer=None):
    r"""
    Print an import statement.

    INPUT:

    - ``module`` -- the name of a module

    - ``name`` -- the name of the object to import

    - ``lazy`` -- a boolean: whether to print a lazy import statement

    EXAMPLES::

        sage: sage.misc.dev_tools.print_import_statement('sage.misc.dev_tools', 'print_import_statement', False)
        from sage.misc.dev_tools import print_import_statement
        sage: sage.misc.dev_tools.print_import_statement('sage.misc.dev_tools', 'print_import_statement', True)
        lazy_import('sage.misc.dev_tools', 'print_import_statement')
    """
    if lazy:
        print_or_update("lazy_import('%s', '%s')"%(module, name), answer)
    else:
        print_or_update("from %s import %s"%(module, name), answer)

def import_statements(*objects, **options):
    """
    Print import statements for the given objects.

    INPUT:

    - ``*objects`` -- a sequence of objects.

    - ``lazy`` -- a boolean (default: ``False``)
      Whether to print a lazy import statement.

    - ``verbose`` -- a boolean (default: ``True``)
      Whether to print information in case of ambiguity.

    - ``answer_as_str`` -- a boolean (default: ``False``)
      If ``True`` return a string instead of printing the statement.

    EXAMPLES::

        sage: import_statements(WeylGroup, lazy_attribute)
        from sage.combinat.root_system.weyl_group import WeylGroup
        from sage.misc.lazy_attribute import lazy_attribute

        sage: import_statements(IntegerRing)
        from sage.rings.integer_ring import IntegerRing

    If ``lazy`` is True, then :func:`lazy_import` statements are
    displayed instead::

        sage: import_statements(WeylGroup, lazy_attribute, lazy=True)
        from sage.misc.lazy_import import lazy_import
        lazy_import('sage.combinat.root_system.weyl_group', 'WeylGroup')
        lazy_import('sage.misc.lazy_attribute', 'lazy_attribute')

    In principle, the function should also work on object which are instances.
    In case of ambiguity, one or two warning lines are printed::

        sage: import_statements(RDF)
        from sage.rings.real_double import RDF

        sage: import_statements(ZZ)
          ** Warning **: several names for that object: Z, ZZ
        from sage.rings.integer_ring import Z

        sage: import_statements(euler_phi)
        from sage.rings.arith import euler_phi

        sage: import_statements(x)
          ** Warning **: several modules for that object: sage.all_cmdline, sage.calculus.predefined
        from sage.calculus.predefined import x

    If you don't like the warning you can disable them with the option ``verbose``::

        sage: import_statements(ZZ, verbose=False)
        from sage.rings.integer_ring import Z

        sage: import_statements(x, verbose=False)
        from sage.calculus.predefined import x

    If the object has several names, an other way to get the import
    statement you expect is to use a string instead of the object::

        sage: import_statements(cached_function)
          ** Warning **: several names for that object: CachedFunction, cached_function
        from sage.misc.cachefunc import CachedFunction

        sage: import_statements('cached_function')
        from sage.misc.cachefunc import cached_function
        sage: import_statements('Z')
        from sage.rings.integer_ring import Z


    Specifying a string is also useful for objects that are not
    imported in the Sage interpreter namespace by default. In this
    case, an object with that name is looked up in all the modules
    that have been imported in this session::

        sage: print_import_statement
        Traceback (most recent call last):
        ...
        NameError: name 'print_import_statement' is not defined

        sage: import_statements("print_import_statement")
        from sage.misc.dev_tools import print_import_statement

    We test different object which have no appropriate answer::

        sage: import_statements('my_tailor_is_rich')
        Traceback (most recent call last):
        ...
        ValueError: no import statement for my_tailor_is_rich
        sage: import_statements(5)
        Traceback (most recent call last):
        ...
        ValueError: no import statement for 5

    We test that it behaves well with lazy imported objects (:trac:`14767`)::

        sage: import_statements(NN)
        from sage.rings.semirings.non_negative_integer_semiring import NN
        sage: import_statements('NN')
        from sage.rings.semirings.non_negative_integer_semiring import NN
    """
    import inspect, sys, re
    import sage.all
    from sage.misc import sageinspect
    from sage.misc.flatten import flatten
    from sage.misc.lazy_import import LazyImport

    answer_as_str = options.get("answer_as_str",False)
    if answer_as_str:
        answer = []
    else:
        answer = None

    lazy = options.get("lazy", False)
    verbose = options.get("verbose", True)
    if lazy:
        print_or_update("from sage.misc.lazy_import import lazy_import", answer)

    for obj in objects:
        # if obj is a string use it has a name and look for an object
        if isinstance(obj, str):
            name = obj
            obj = sage.all.__dict__.get(name)
            if obj is None:
                # Look for the first module which contains that name.
                # TODO: improve this heuristic.
                for module in sys.modules.values():
                    if hasattr(module, '__dict__') and name in module.__dict__:
                        obj = module.__dict__[name]
                        break
                else:
                    raise ValueError("no import statement for %s"%name)

        else:
            name = None


        if isinstance(obj, LazyImport):
            obj = obj._get_object()

        # Case 1: the object is a module
        if inspect.ismodule(obj):
            if lazy:
                print_or_update("lazy_import('%s')"%obj.__name__, answer)
            else:
                print_or_update("import %s"%obj.__name__, answer)
            continue

        # Case 2: the object is defined in its module
        module = None
        if sageinspect.isclassinstance(obj):
            module = obj.__class__.__module__
        elif hasattr(obj, '__module__') and obj.__module__:
            module = obj.__module__

        if module:
            d = sys.modules[module].__dict__
            names = None
            if name is None:
                names = sorted(key for key in d if d[key] is obj)
            elif name in d:
                names = [name]
            if names:
                if verbose and len(names) > 1:
                    print "  ** Warning **: several names for that object: %s"%', '.join(names)
                print_import_statement(module, names[0], lazy, answer)
                continue


        # Case 3: search for this object in all modules
        names = {} # dictionnary: module -> names of the object in that module
        for module in sys.modules:
            if module != '__main__' and hasattr(sys.modules[module],'__dict__'):
                d = sys.modules[module].__dict__

                if name is not None:
                    if name in d and d[name] is obj:
                        names[module] = [name]
                else:
                    n = [key for key in d if d[key] is obj]
                    if n:
                        names[module] = n

        all_names = sorted(set(flatten(names.values())))
        if len(all_names) == 0:
            raise ValueError("no import statement for %s"%obj)
        elif verbose and len(all_names) > 1:
            print "  ** Warning **: several names for that object:",
            print ", ".join(sorted(all_names))

        modules = sorted(flatten(names),cmp=module_names_cmp)
        if verbose and len(modules) > 1:
            print "  ** Warning **: several modules for that object:",
            print ", ".join(modules[:4]),
            if len(modules) > 4:
                print "..."
            else:
                print

        # Case 4: if the object is a class instance, we look for a
        # module where it is instanciated
        if sageinspect.isclassinstance(obj):
            names_pattern = dict((name,re.compile("^%s\ *="%name, re.MULTILINE)) for name in all_names)

            # if obj has a module we try to put the .all of that object at the beginig
            # (it at least solves the problem for CC and CIF)
            if hasattr(obj, '__module__'):
                module = obj.__module__
                try:
                    i = module.rindex('.')
                    new_module = module[:i] + '.all'
                except ValueError:
                    new_module = module

                if new_module in modules:
                    i = modules.index(new_module)
                    del modules[i]
                    modules.insert(0,new_module)

            for module in modules:
                sources = sageinspect.sage_getsource(sys.modules[module])
                for name in names[module]:
                    if names_pattern[name].search(sources):
                        break
                else:
                    continue
                break
        else:
            module = modules[0]
            name = names[module][0]

        if name is not None:
            print_import_statement(module, name, lazy, answer)
        else:
            raise ValueError("no import statement for %s"%obj)

    if answer is not None:
        return '\n'.join(answer)

def which_import_statements_fail():
    r"""
    Run import statements on *all* objects and print the one for which the
    import fails or return "from sage.all import my_object".

    The returne value is a couple of lists. The first one corresponds to the
    import_statements that answer "from sage.all import XXX" and the second
    corresponds to the list of import_statements that fail.

    TESTS::

        sage: from sage.misc.dev_tools import which_import_statements_fail
        sage: sage_all, errors = which_import_statements_fail() # long time
        ** Warning **: several modules for that object: sage.all, sage.all_cmdline
        ...
        sage: errors          # long time
        ['maxima_calculus']
        sage: len(sage_all)   # long time
        62
    """
    import sage.all
    errors = []
    sage_all = []
    for name in sage.all.__dict__.iterkeys():
        try:
            to_test = import_statements(name, answer_as_str=True)
        except ValueError:
            errors.append(name)
            continue

        if to_test.startswith('from sage.all import'):
            sage_all.append(name)
            continue
        try:
            exec import_statements(name, answer_as_str=True)
        except ImportError:
            errors.append(name)
    return sage_all, errors
