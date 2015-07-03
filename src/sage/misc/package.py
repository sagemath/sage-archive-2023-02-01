r"""
Sage package management commands

A Sage package has the extension .spkg. It is a tarball that is
(usually) bzip2 compressed that contains arbitrary data and an
spkg-install file. An Sage package typically has the following
components:


-  spkg-install - shell script that is run to install the package

-  Sage.txt - file that describes how the package was made, who
   maintains it, etc.

-  sage - directory with extra patched version of files that needed
   during the install


Use the :func:`install_package` command to install a new package, and use
:func:`optional_packages` to list all optional packages available. The
:func:`upgrade` command upgrades all *standard* packages - there is no
auto-upgrade command for optional packages.

All package management can also be done via the Sage command line.

Packages available
------------------

**Standard packages:**

{STANDARD_PACKAGES}

**Optional packages:**

{OPTIONAL_PACKAGES}

**Experimental packages:**

{EXPERIMENTAL_PACKAGES}

Functions
---------
"""

def _list_to_table(list_of_packages):
    r"""
    Helper function returning a ReST table from a list of strings.

    The entries are sorted vertically.

    INPUT:

    - ``list_of_packages`` -- a list

    EXAMPLE::

        sage: print sage.misc.package._list_to_table([str(x) for x in range(10)])
        .. csv-table::
            :class: contentstable
            :widths: 20, 20, 20, 20, 20
            :delim: |
        <BLANKLINE>
                ``0`` | ``2`` | ``4`` | ``6`` | ``8``
                ``1`` | ``3`` | ``5`` | ``7`` | ``9``
        <BLANKLINE>

    Check that the local list of packages matches the online list. Standard
    packages::

        sage: from sage.misc.package import _STANDARD_PACKAGES, standard_packages
        sage: a,b = standard_packages() # optional internet
        sage: set(a+b).symmetric_difference(_STANDARD_PACKAGES) # optional internet
        set()

    Optional packages::

        sage: from sage.misc.package import _OPTIONAL_PACKAGES, optional_packages
        sage: a,b = optional_packages() # optional internet
        sage: set(a+b).symmetric_difference(_OPTIONAL_PACKAGES) # optional internet
        set()

    Experimental packages::

        sage: from sage.misc.package import _EXPERIMENTAL_PACKAGES, experimental_packages
        sage: a,b = experimental_packages() # optional internet
        sage: set(a+b).symmetric_difference(_EXPERIMENTAL_PACKAGES) # optional internet
        set()
    """
    from string import join
    s = (".. csv-table::\n"
                "    :class: contentstable\n"
                "    :widths: 20, 20, 20, 20, 20\n"
                "    :delim: |\n\n")
    length = len(list_of_packages)
    width = 5
    height = (length+width-1)//width

    list_of_packages = sorted(["``"+p+"``" if p else p
                               for p in list_of_packages])

    list_of_packages.sort()
    list_of_packages.extend(['']*width)
    for l in range(height):
        s += "        "+join(list_of_packages[l::height][:width], ' | ')+"\n"

    return s

_STANDARD_PACKAGES = ['atlas', 'backports_ssl_match_hostname', 'boehm_gc',
                      'boost_cropped', 'bzip2', 'cddlib', 'cephes', 'certifi', 'cliquer',
                      'combinatorial_designs', 'conway_polynomials', 'cvxopt', 'cython', 'dateutil',
                      'docutils', 'ecl', 'eclib', 'ecm', 'elliptic_curves', 'fflas_ffpack', 'flint',
                      'flintqs', 'freetype', 'gap', 'gd', 'gdmodule', 'genus2reduction', 'gf2x',
                      'gfan', 'git', 'givaro', 'glpk', 'graphs', 'gsl', 'iconv', 'iml', 'ipython',
                      'jinja2', 'jmol', 'jsonschema', 'lcalc', 'libfplll', 'libgap', 'libgd',
                      'libpng', 'linbox', 'lrcalc', 'm4ri', 'm4rie', 'markupsafe', 'mathjax',
                      'matplotlib', 'maxima', 'mercurial', 'mistune', 'mpc', 'mpfi', 'mpfr', 'mpmath',
                      'ncurses', 'networkx', 'ntl', 'numpy', 'palp', 'pari', 'pari_galdata',
                      'pari_seadata_small', 'patch', 'pexpect', 'pil', 'pillow', 'pip', 'pkgconf',
                      'pkgconfig', 'planarity', 'polybori', 'polytopes_db', 'ppl', 'pycrypto',
                      'pygments', 'pynac', 'pyparsing', 'python', 'pyzmq', 'r', 'ratpoints',
                      'readline', 'rpy2', 'rubiks', 'rw', 'sage', 'sage_root', 'sage_scripts',
                      'sagenb', 'sagetex', 'scipy', 'setuptools', 'singular', 'six', 'sphinx',
                      'sqlalchemy', 'sqlite', 'symmetrica', 'sympow', 'sympy', 'tachyon', 'tornado',
                      'zeromq', 'zlib', 'zn_poly']

_OPTIONAL_PACKAGES = ['PyQt_x11', 'TOPCOM', 'arb', 'beautifulsoup', 'benzene',
                      'biopython', 'bliss', 'brian', 'buckygen', 'cbc', 'ccache', 'chomp',
                      'cluster_seed', 'coxeter3', 'cryptominisat', 'cunningham_tables', 'd3js',
                      'database_cremona_ellcurve', 'database_gap', 'database_jones_numfield',
                      'database_kohel', 'database_odlyzko_zeta', 'database_pari',
                      'database_stein_watkins', 'database_stein_watkins_mini',
                      'database_symbolic_data', 'dot2tex', 'extra_docs', 'gambit', 'gap_packages',
                      'gcc', 'gdb', 'giac', 'giacpy', 'ginv', 'git_trac', 'gmp', 'gnuplotpy', 'guppy',
                      'java3d', 'kash3', 'knoboo', 'libogg', 'libtheora', 'lidia', 'lie', 'lrslib',
                      'mcqd', 'modular_decomposition', 'mpi4py', 'mpir', 'nauty', 'normaliz', 'nose',
                      'nzmath', 'openmpi', 'openssl', 'ore_algebra', 'p_group_cohomology', 'phc',
                      'plantri', 'pybtex', 'python2', 'python3', 'pyx', 'qhull', 'sage_mode', 'scons',
                      'sip', 'termcap', 'threejs', 'tides', 'topcom', 'trac']


_EXPERIMENTAL_PACKAGES = ['4ti2', 'PyQt4', 'PyVTK', 'QScintilla2', 'asymptote',
                          'autotools', 'bison', 'cadabra', 'clapack', 'clisp', 'cmake', 'compilerwrapper',
                          'csdp', 'dvipng', 'ets', 'fes', 'flex', 'fricas', 'gnofract4d', 'gnuplot',
                          'graphviz', 'latte_int', 'libcprops', 'libjpeg', 'libsigsegv', 'macaulay2',
                          'mayavi', 'meataxe', 'modglue', 'mpich2', 'numarray', 'numeric', 'openopt',
                          'pcre', 'phcpack', 'polymake', 'processing', 'pygame', 'pygsl', 'pygtk',
                          'pynifti', 'pyqt', 'pyrexembed', 'qasm', 'qepcad', 'quantlib', 'quantlib_swig',
                          'reallib3_linux', 'sandpile', 'scitools++', 'semigroupe', 'simpqs', 'sip',
                          'soya', 'soya_cvs', 'superlu', 'surf', 'valgrind', 'vtk_meta', 'wxPython',
                          'yafray', 'yassl']


__doc__ = __doc__.format(STANDARD_PACKAGES     =_list_to_table(_STANDARD_PACKAGES),
                         OPTIONAL_PACKAGES     =_list_to_table(_OPTIONAL_PACKAGES),
                         EXPERIMENTAL_PACKAGES =_list_to_table(_EXPERIMENTAL_PACKAGES))

import os

__installed_packages = None

def install_all_optional_packages(force=True, dry_run=False):
    r"""
    Install all available optional spkg's in the official Sage spkg
    repository.  Returns a list of all spkg's that *fail* to install.

    INPUT:

    - ``force`` -- bool (default: ``True``); whether to force
      reinstall of spkg's that are already installed.

    - ``dry_run`` -- bool (default: ``False``); if ``True``, just list
      the packages that would be installed in order, but don't
      actually install them.

    OUTPUT:

        list of strings

    .. NOTE::

        This is designed mainly for testing purposes.  This also
        doesn't do anything with respect to dependencies -- the
        packages are installed in alphabetical order.  Dependency
        issues will be dealt with in a future version.


    AUTHOR:

        -- William Stein (2008-12)

    EXAMPLES::

        sage: sage.misc.package.install_all_optional_packages(dry_run=True)  # optional - internet
        Installing ...
        []
    """
    # Get list of all packages from the package server
    installed, not_installed = optional_packages()
    if force:
        v = installed + not_installed
    else:
        v = not_installed
    failed = []

    # sort the list of packages in alphabetical order
    v.sort()

    # install each one
    for pkg in v:
        try:
            print("Installing {}...".format(pkg))
            if not dry_run:
                # only actually do the install of the package if dry_run is False
                install_package(pkg, force=force)
        except ValueError as msg:
            # An error occurred -- catch exception and record this in
            # the list of failed installed.
            print("*"*70)
            print("FAILED to install '{}'".format(pkg))
            print("*"*70)
            failed.append(pkg)
    return failed

def install_package(package=None, force=False):
    """
    Install a package or return a list of all packages that have been
    installed into this Sage install.

    You must have an internet connection. Also, you will have to
    restart Sage for the changes to take affect.

    It is not needed to provide the version number.

    INPUT:

    -  ``package`` -- optional; if specified, install the
       given package. If not, list all installed packages.

    - ``force`` -- boolean (default: ``False``); if ``True``, reinstall
      the package given if it is already installed. (Otherwise, an already
      installed package doesn't get reinstalled, as with 'sage -i ...'.)

    EXAMPLES:

    With no arguments, list the installed packages::

        sage: install_package()
        [...'atlas...'python...]

    With an argument, install the named package::

        sage: install_package('chomp')  # not tested
        Attempting to download package chomp-20100213.p2
        ...

    IMPLEMENTATION:

    Calls 'sage -f ...' to (re)install the package if a package name is
    given. If no package name is given, simply list the contents of
    ``spkg/installed``.

    .. seealso:: :func:`optional_packages`, :func:`upgrade`
    """
    global __installed_packages
    if os.uname()[0][:6] == 'CYGWIN' and package is not None:
        print("install_package may not work correctly under Microsoft Windows")
        print("since you can't change an opened file.  Quit all")
        print("instances of Sage and use 'sage -i {}' instead or".format(package))
        print("use the force option to install_package().")

    if package is None:
        if __installed_packages is None:
            import sage.env
            __installed_packages = sorted(os.listdir(sage.env.SAGE_SPKG_INST))
        return __installed_packages



    import stopgap
    stopgap.stopgap("The Sage function 'install_packages' is currently broken: it does not correctly install new packages. Please use 'sage -i {}' from a shell prompt instead.".format(package), 16759)
    return



    # Get full package name / list of candidates:
    if force:
        # Also search packages already installed.
        S = [P for P in standard_packages()[0] if P.startswith(package)]
        O = [P for P in optional_packages()[0] if P.startswith(package)]
        E = [P for P in experimental_packages()[0] if P.startswith(package)]
    else:
        S,O,E = [], [], []
    S.extend([P for P in standard_packages()[1] if P.startswith(package)])
    O.extend([P for P in optional_packages()[1] if P.startswith(package)])
    E.extend([P for P in experimental_packages()[1] if P.startswith(package)])
    L = S+O+E
    if len(L) > 1:
        if force:
            print("Possible package names starting with '{}' are:".format(package))
        else:
            print("Possible names of non-installed packages starting with '{}':".format(package))
        for P in L:
            print(" ", P)
        raise ValueError("There is more than one package name starting with '{}'. Please specify!".format(package))
    if len(L) == 0:
        if not force:
            if is_package_installed(package):
                raise ValueError("Package is already installed. Try install_package('{}',force=True)".format(package))
        raise ValueError("There is no package name starting with '{}'.".format(package))
    # len(L) == 1, i.e. exactly one package matches the given one.
    os.system('sage -f "{}"'.format(L[0]))
    __installed_packages = None


def is_package_installed(package):
    """
    Return true if ``package`` is installed.

    EXAMPLES::

        sage: is_package_installed('pari')
        True

    Giving just the beginning of the package name is not good enough::

        sage: is_package_installed('matplotli')
        False

    Otherwise, installing "pillow" will cause this function to think
    that "pil" is installed, for example.
    """
    return any(p.split('-')[0] == package for p in install_package())

def package_versions(package_type, local=False):
    r"""
    Return version information for each Sage package.

    INPUT:

    - ``package_type`` (string) -- one of `"standard"`, `"optional"` or
      `"experimental"`

    - ``local`` (boolean) -- only query local data (no internet needed)

    For packages of the given type, return a dictionary whose entries
    are of the form ``'package': (installed, latest)``, where
    ``installed`` is the installed version (or ``None`` if not
    installed) and ``latest`` is the latest available version. If the
    package has a directory in ``SAGE_ROOT/build/pkgs/``, then
    ``latest`` is determined by the file ``package-version.txt`` in
    that directory.  If ``local`` is False, then Sage's servers are
    queried for package information.

    EXAMPLES::

        sage: std = package_versions('standard', local=True)
        sage: 'gap' in std
        True
        sage: std['zn_poly']
        ('0.9.p11', '0.9.p11')
    """
    if package_type not in ['standard','optional','experimental']:
        raise ValueError("'package_type' must be one of 'standard','optional','experimental'.")

    cmd = 'sage-list-packages {} --dump'.format(package_type)
    if local:
        cmd += " --local"
    X = os.popen(cmd).read().split('\n')[:-1]

    versions = {}
    for line in X:
        line = line.split(' ')
        installed = line[2]
        if installed == 'not_installed': installed = None
        versions[line[0]] = (installed, line[1])
    return versions

def _package_lists_from_sage_output(package_type, local=False):
    r"""
    Helper function for :func:`standard_packages`, :func:`optional_packages` and
    :func:`experimental_packages`.

    INPUT:

    - ``package_type`` (string) -- one of `"standard"`, `"optional"` or
      `"experimental"`

    - ``local`` (boolean) -- only query local data (no internet needed)

    OUTPUT:

    The function returns a pair of lists ``(installed,not_installed)``
    with the corresponding packages' name, sorted alphabetically. If
    ``local`` is ``True``, then the list of all packages is downloaded
    from the server; otherwise, the list is extracted from the
    packages in ``SAGE_ROOT/build/pkgs/`

    EXAMPLE::

        sage: from sage.misc.package import standard_packages
        sage: installed, not_installed = standard_packages()  # optional internet

    Local check that all standard packages are installed::

        sage: from sage.misc.package import _package_lists_from_sage_output
        sage: installed, not_installed = _package_lists_from_sage_output('standard',local=True)
        sage: bool(not_installed)
        False
        sage: 'glpk' in installed
        True
    """
    installed     = []
    not_installed = []
    versions = package_versions(package_type, local)
    for p in versions:
        if versions[p][0] is None:
            not_installed.append(p)
        else:
            installed.append(p)

    return sorted(installed), sorted(not_installed)

def standard_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed standard packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed standard packages (as a list)

    -  NOT installed standard packages (as a list)

    Use ``install_package(package_name)`` to install or
    re-install a given package.

    .. seealso:: :func:`install_package`, :func:`upgrade`

    EXAMPLE::

        sage: from sage.misc.package import standard_packages
        sage: installed, not_installed = standard_packages() # optional internet
        sage: installed[0], installed[-1]                    # optional internet
        ('atlas', 'zn_poly')
        sage: 'mercurial' in not_installed                   # optional internet
        True

    """
    return _package_lists_from_sage_output('standard')

def optional_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed optional packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed optional packages (as a list)

    -  NOT installed optional packages (as a list)


    Use ``install_package(package_name)`` to install or
    re-install a given package.

    .. seealso:: :func:`install_package`, :func:`upgrade`

    EXAMPLE::

        sage: from sage.misc.package import optional_packages
        sage: installed, not_installed = optional_packages() # optional internet
        sage: min(installed+not_installed)                   # optional internet
        'PyQt_x11'
        sage: max(installed+not_installed)                   # optional internet
        'trac'
    """
    return _package_lists_from_sage_output('optional')

def experimental_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed experimental packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed experimental packages (as a list)

    -  NOT installed experimental packages (as a list)


    Use ``install_package(package_name)`` to install or
    re-install a given package.

    .. seealso:: :func:`install_package`, :func:`upgrade`

    EXAMPLE::

        sage: from sage.misc.package import experimental_packages
        sage: installed, not_installed = experimental_packages() # optional internet
        sage: min(installed+not_installed)                   # optional internet
        '4ti2'
        sage: max(installed+not_installed)                   # optional internet
        'yassl'
    """
    return _package_lists_from_sage_output('experimental')

#################################################################
# Upgrade to latest version of Sage
#################################################################


def upgrade():
    """
    Download and build the latest version of Sage.

    You must have an internet connection. Also, you will have to
    restart Sage for the changes to take affect.

    This upgrades to the latest version of core packages (optional
    packages are not automatically upgraded).

    This will not work on systems that don't have a C compiler.

    .. seealso:: :func:`install_package`, :func:`optional_packages`
    """
    global __installed_packages
    if os.uname()[0][:6] == 'CYGWIN':
        print("Upgrade may not work correctly under Microsoft Windows")
        print("since you can't change an opened file.  Quit all")
        print("instances of Sage and use 'sage -upgrade' instead.")
        return

    os.system('sage -upgrade')
    __installed_packages = None
    print("You should quit and restart Sage now.")


def package_mesg(package_name):
    mesg  = 'To install the package %s type install_package("%s")\n'%(package_name, package_name)
    mesg += 'at the sage prompt.  Note, the version number might\n'
    mesg += 'change; if so, type optional_packages() to see a list \n'
    mesg += 'of possibilities.   All this requires an internet connection.'
    mesg += 'For more help, type optional_packages?'
    return mesg

