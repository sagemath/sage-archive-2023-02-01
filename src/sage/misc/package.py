r"""
Listing Sage packages

This module can be used to see which Sage packages are installed
and which packages are available for installation.

For more information about creating Sage packages, see
the "Packaging Third-Party Code" section of the
Sage Developer's Guide.

Actually installing the packages should be done via the command
line, using the following commands:

- ``sage -i PACKAGE_NAME`` -- install the given package

- ``sage -f PACKAGE_NAME`` -- re-install the given package, even if it
  was already installed
  
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

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

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

def install_package(package=None, force=None):
    """
    This function is obsolete. Run ``sage -i PKGNAME`` from a shell
    to install a package. Use the function :func:`installed_packages`
    to list all installed packages.

    TESTS::

        sage: install_package()
        doctest:...: DeprecationWarning: use installed_packages() to list all installed packages
        See http://trac.sagemath.org/16759 for details.
        [...'atlas...'python...]
        sage: install_package("autotools")
        Traceback (most recent call last):
        ...
        NotImplementedError: installing Sage packages using 'install_package()' is obsolete.
        Run 'sage -i autotools' from a shell prompt instead
    """
    if package is not None:
        # deprecation(16759, ...)
        raise NotImplementedError("installing Sage packages using 'install_package()' is obsolete.\nRun 'sage -i {}' from a shell prompt instead".format(package))

    from sage.misc.superseded import deprecation
    deprecation(16759, "use installed_packages() to list all installed packages")
    return installed_packages()


def installed_packages():
    """
    Return a list of all installed packages, with version numbers.

    EXAMPLES::

        sage: installed_packages()
        [...'atlas...'python...]

    .. seealso::

        :func:`standard_packages`, :func:`optional_packages`, :func:`experimental_packages`
    """
    from sage.env import SAGE_SPKG_INST
    return sorted(os.listdir(SAGE_SPKG_INST))


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
    return any(p.split('-')[0] == package for p in installed_packages())

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

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

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

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

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

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    EXAMPLE::

        sage: from sage.misc.package import experimental_packages
        sage: installed, not_installed = experimental_packages() # optional internet
        sage: min(installed+not_installed)                   # optional internet
        '4ti2'
        sage: max(installed+not_installed)                   # optional internet
        'yassl'
    """
    return _package_lists_from_sage_output('experimental')

def upgrade():
    """
    Obsolete function, run 'sage --upgrade' from a shell prompt instead.

    TESTS::

        sage: upgrade()
        Traceback (most recent call last):
        ...
        NotImplementedError: upgrading Sage using 'upgrade()' is obsolete.
        Run 'sage --upgrade' from a shell prompt instead
    """
    # deprecation(16759, ..)
    raise NotImplementedError("upgrading Sage using 'upgrade()' is obsolete.\nRun 'sage --upgrade' from a shell prompt instead")

class PackageNotFoundError(RuntimeError):
    """
    This class defines the exception that should be raised when a
    function, method, or class cannot detect a Sage package that it
    depends on.

    This exception should be raised with a single argument, namely
    the name of the package.

    When an ``PackageNotFoundError`` is raised, this means one of the
    following:

    - The required optional package is not installed.

    - The required optional package is installed, but the relevant
      interface to that package is unable to detect the package.

    EXAMPLES::

        sage: from sage.misc.package import PackageNotFoundError
        sage: raise PackageNotFoundError("my_package")
        Traceback (most recent call last):
        ...
        PackageNotFoundError: the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell
    """
    def __str__(self):
        """
        Return the actual error message.

        EXAMPLES::

            sage: from sage.misc.package import PackageNotFoundError
            sage: str(PackageNotFoundError("my_package"))
            "the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell"
        """
        return ("the package {0!r} was not found. "
            "You can install it by running 'sage -i {0}' in a shell"
            .format(self.args[0]))
