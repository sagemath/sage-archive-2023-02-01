# -*- coding: utf-8 -*-
r"""
Testing for features of the environment at runtime

A computation can require a certain package to be installed in the runtime
environment. Abstractly such a package describes a :class`Feature` which can
be tested for at runtime. It can be of various kinds, most prominently an
:class:`Executable` in the PATH or an additional package for some installed
system such as a :class:`GapPackage`.

AUTHORS:

- Julian Rüth (2016-04-07): Initial version

EXAMPLES:

Some generic features are available for common cases. For example, to
test for the existence of a binary, one can use an :class:`Executable`
feature::

    sage: from sage.misc.feature_test import GapPackage, Executable
    sage: Executable(name="sh", executable="sh").is_present()
    FeatureTestResult('sh', True)

Here we test whether the grape GAP package is available::

    sage: GapPackage("grape", spkg="gap_packages").is_present() # optional: gap_packages
    FeatureTestResult('GAP package grape', True)

Note that a :class:`FeatureTestResult` acts like a bool in most contexts::

    sage: if Executable(name="sh", executable="sh").is_present(): "present."
    'present.'

When one wants to raise an error if the feature is not available, one
can use the ``require`` method::

    sage: Executable(name="sh", executable="sh").require()

    sage: Executable(name="random", executable="randomOochoz6x", spkg="random", url="http://rand.om").require()
    Traceback (most recent call last):
    ...
    FeatureNotPresentError: random is not available.
    No `randomOochoz6x` found on PATH.
    To install random you can try to run `sage -i random`.
    Further installation instructions might be available at http://rand.om.

As can be seen above, features try to produce helpful error messages.
"""
from sage.misc.cachefunc import cached_method

class Feature(object):
    r"""
    A feature of the runtime environment

    Overwrite :meth:`is_present` to add feature checks.

    EXAMPLES::

        sage: from sage.misc.feature_test import GapPackage
        sage: GapPackage("grape", spkg="gap_packages") # indirect doctest
        Feature('GAP package grape')
    """
    def __init__(self, name, spkg = None, url = None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import GapPackage, Feature
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), Feature) # indirect doctest
            True
        """
        self.name = name;
        self.spkg = spkg
        self.url = url

    def is_present(self, explain=False):
        r"""
        Return whether the feature is present.

        OUTPUT:

        A :class:`FeatureTestResult` which can be used as a boolean and
        contains additional information about the feature test.

        EXAMPLES::

            sage: from sage.misc.feature_test import GapPackage
            sage: GapPackage("grape", spkg="gap_packages").is_present() # optional: gap_packages
            FeatureTestResult('GAP package grape', True)
            sage: GapPackage("NOT_A_PACKAGE", spkg="gap_packages").is_present()
            FeatureTestResult('GAP package NOT_A_PACKAGE', False)
        """
        return FeatureTestResult(self, True)

    def require(self):
        r"""
        Raise a :class:`FeatureNotPresentError` if the feature is not present.

        EXAMPLES::

            sage: from sage.misc.feature_test import GapPackage
            sage: GapPackage("ve1EeThu").require()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package ve1EeThu is not available.
            `TestPackageAvailability("ve1EeThu")` evaluated to `fail` in GAP.
        """
        presence = self.is_present()
        if not presence:
            raise FeatureNotPresentError(self, presence.reason, presence.resolution)

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        EXAMPLES::

            sage: from sage.misc.feature_test import GapPackage
            sage: GapPackage("grape") # indirect doctest
            Feature('GAP package grape')
        """
        return 'Feature({name!r})'.format(name=self.name)

    def resolution(self):
        r"""
        Return a suggestion on how to make :meth:`is_present` pass if it did not
        pass.

        EXAMPLES::

            sage: from sage.misc.feature_test import Executable
            sage: Executable(name="CSDP", spkg="csdp", executable="theta", url="http://github.org/dimpase/csdp").resolution()
            'To install CSDP you can try to run `sage -i csdp`.\nFurther installation instructions might be available at http://github.org/dimpase/csdp.'
        """
        return "\n".join(filter(None,[
            "To install {feature} you can try to run `sage -i {spkg}`.".format(feature=self.name, spkg=self.spkg) if self.spkg else "",
            "Further installation instructions might be available at {url}.".format(url=self.url) if self.url else ""])) or None

class FeatureNotPresentError(RuntimeError):
    r"""
    A missing feature error.

    EXAMPLES::

        sage: from sage.misc.feature_test import Feature, FeatureTestResult
        sage: class Missing(Feature):
        ....:     def is_present(self): return FeatureTestResult(self, False)

        sage: Missing(name="missing").require()
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: missing is not available.
    """
    def __init__(self, feature, reason=None, resolution=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Feature, FeatureNotPresentError, FeatureTestResult
            sage: class Missing(Feature):
            ....:     def is_present(self): return FeatureTestResult(self, False)

            sage: try:
            ....:     Missing(name="missing").require() # indirect doctest
            ....: except FeatureNotPresentError: pass
        """
        self.feature = feature
        self.reason = reason
        self.resolution = resolution

    def __str__(self):
        r"""
        Return the error message.

        EXAMPLES::

            sage: from sage.misc.feature_test import GapPackage
            sage: GapPackage("gapZuHoh8Uu").require()   # indirect doctest
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package gapZuHoh8Uu is not available.
            `TestPackageAvailability("gapZuHoh8Uu")` evaluated to `fail` in GAP.
        """
        return "\n".join(filter(None,(
            "{feature} is not available.".format(feature=self.feature.name),
            self.reason,
            self.resolution
            )))
        return

class FeatureTestResult(object):
    r"""
    The result of a :method:`Feature.is_present` call.

    Behaves like a boolean with some extra data which may explain why a feature
    is not present and how this may be resolved.

    EXAMPLES::

        sage: from sage.misc.feature_test import GapPackage
        sage: presence = GapPackage("NOT_A_PACKAGE").is_present(); presence # indirect doctest
        FeatureTestResult('GAP package NOT_A_PACKAGE', False)
        sage: bool(presence)
        False

    Explanatory messages might be available as ``reason`` and
    ``resolution``::

        sage: presence.reason
        '`TestPackageAvailability("NOT_A_PACKAGE")` evaluated to `fail` in GAP.'
        sage: presence.resolution is None
        True

    If a feature is not present, ``resolution`` defaults to
    ``feature.resolution()`` if this is defined. If you do not want to use this
    default you need explicitly set ``resolution`` to a string::

        sage: from sage.misc.feature_test import FeatureTestResult
        sage: package = GapPackage("NOT_A_PACKAGE", spkg="no_package")
        sage: FeatureTestResult(package, True).resolution

        sage; FeatureTestResult(package, False).resolution

        sage: FeatureTestResult(package, False, resolution="rtm").resolution
        'rtm'
    """
    def __init__(self, feature, is_present, reason=None, resolution=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Executable, FeatureTestResult
            sage: isinstance(Executable(name="sh", executable="sh").is_present(), FeatureTestResult)
            True
        """
        self.feature = feature
        self.is_present = is_present
        self.reason = reason
        self.resolution = resolution
        if not is_present and self.resolution is None and hasattr(feature, "resolution"):
            self.resolution = self.feature.resolution()

    def __nonzero__(self):
        r"""
        Whether the tested :class:`Feature` is present.

        TESTS::

            sage: from sage.misc.feature_test import Feature, FeatureTestResult
            sage: bool(FeatureTestResult(Feature("SomePresentFeature"), True)) # indirect doctest
            True
            sage: bool(FeatureTestResult(Feature("SomeMissingFeature"), False))
            False
        """
        return bool(self.is_present)

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Feature, FeatureTestResult
            sage: FeatureTestResult(Feature("SomePresentFeature"), True) # indirect doctest
            FeatureTestResult('SomePresentFeature', True)
        """
        return "FeatureTestResult({feature!r}, {is_present!r})".format(feature=self.feature.name, is_present=self.is_present)

class Executable(Feature):
    r"""
    A feature describing an executable in the PATH.

    .. NOTE::

        Overwrite :meth:`is_functional` if you also want to check whether
        the executable shows proper behaviour.

        Calls to :meth:`is_present` are cached. You might want to cache the
        :class:`Executable` object to prevent unnecessary calls to the
        executable.

    EXAMPLES::

        sage: from sage.misc.feature_test import Executable
        sage: Executable(name="sh", executable="sh").is_present()
        FeatureTestResult('sh', True)
    """
    def __init__(self, name, executable, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Executable
            sage: isinstance(Executable(name="sh", executable="sh"), Executable)
            True
        """
        Feature.__init__(self, name=name, spkg=spkg, url=url)
        self.executable = executable

    @cached_method
    def is_present(self):
        r"""
        Test whether the executable is on the current PATH and functional.

        .. SEEALSO:: :meth:`is_functional`

        EXAMPLES::

            sage: from sage.misc.feature_test import Executable
            sage: Executable(name="sh", executable="sh").is_present()
            FeatureTestResult('sh', True)
        """
        from distutils.spawn import find_executable
        if find_executable(self.executable) is None:
            return FeatureTestResult(self, False, "No `{executable}` found on PATH.".format(executable=self.executable))
        return self.is_functional()

    def is_functional(self):
        r"""
        Return whether an executable in the path is functional.

        EXAMPLES:

        Returns ``True`` unless explicitly overwritten::

            sage: from sage.misc.feature_test import Executable
            sage: Executable(name="sh", executable="sh").is_functional()
            FeatureTestResult('sh', True)
        """
        return FeatureTestResult(self, True)

class StaticFile(Feature):
    r"""
    A :class:`Feature` which describes the presence of a certain file such as a
    database.

    EXAMPLES::

        sage: from sage.misc.feature_test import StaticFile
        sage: StaticFile(name="no_such_file", filename="KaT1aihu", search_path=["/"], spkg="some_spkg", url="http://rand.om").require()
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: no_such_file is not available.
        `KaT1aihu` not found  in any of ['/']
        To install no_such_file you can try to run `sage -i some_spkg`.
        Further installation instructions might be available at http://rand.om.
    """
    def __init__(self, name, filename, search_path, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import StaticFile
            sage: StaticFile(name="null", filename="null", search_path=["/dev"])
            Feature('null')
        """
        Feature.__init__(self, name=name, spkg=spkg, url=url)
        self.filename = filename
        self.search_path = search_path

    def is_present(self):
        r"""
        Whether the static file is present.

           sage: from sage.misc.feature_test import StaticFile
           sage: StaticFile(name="no_such_file", filename="KaT1aihu", search_path=["/"], spkg="some_spkg", url="http://rand.om").is_present()
           FeatureTestResult('no_such_file', False)
        """
        try:
            abspath = self.absolute_path()
            return FeatureTestResult(self, True, reason="Found at `{abspath}`.".format(abspath=abspath))
        except FeatureNotPresentError as e:
            return FeatureTestResult(self, False, reason=e.reason, resolution=e.resolution)

    def absolute_path(self):
        r"""
        The absolute path of the file.

        EXAMPLES::

            sage: from sage.misc.feature_test import DatabaseCremona
            sage: DatabaseCremona().absolute_path() # optional: database_cremona_ellcurve
            '.../local/share/cremona/cremona.db'

        A ``FeatureNotPresentError`` is raised if the file can not be found::

            sage: from sage.misc.feature_test import StaticFile
            sage: StaticFile(name="no_such_file", filename="KaT1aihu", search_path=["/"], spkg="some_spkg", url="http://rand.om").absolute_path()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: no_such_file is not available.
            `KaT1aihu` not found  in any of ['/']
            To install no_such_file you can try to run `sage -i some_spkg`.
            Further installation instructions might be available at http://rand.om.
        """
        for directory in self.search_path:
            import os.path
            import shutil
            path = os.path.join(directory, self.filename)
            if os.path.isfile(path):
                return os.path.abspath(path)
        raise FeatureNotPresentError(self, reason="`{filename}` not found  in any of {search_path}".format(filename=self.filename, search_path=self.search_path), resolution=self.resolution())

class SharedLibrary(Feature):
    r"""
    A :class:`Feature` which describes the presence of a shared library to link
    against from Cython.

    To test the presence of ``name``, the cython compiler is run on
    ``test_code`` and the resulting module is imported.

    EXAMPLES::

        sage: from sage.misc.feature_test import SharedLibrary
        sage: libm_test_code = '''
        ....: #clib m
        ....: cdef extern from "<math.h>":
        ....:     double fabs(double x)
        ....:
        ....: sig_on()
        ....: if fabs(-1) != 1:
        ....:     raise ImportError("libm did not find the correct absolute value for -1")
        ....: sig_off()'''
        sage: libm = SharedLibrary("libm", test_code=libm_test_code, spkg="gcc", url="https://gnu.org")
        sage: libm.is_present()
        FeatureTestResult('libm', True)
    """
    def __init__(self, name, test_code, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import SharedLibrary, LibFES
            sage: any([isinstance(dep, SharedLibrary) for dep in LibFES().dependencies]) # indirect doctest
            True
        """
        Feature.__init__(self, name, spkg=spkg, url=url)
        self.test_code = test_code

    def is_present(self):
        r"""
        Run test code to determine whether the shared library is present.

        EXAMPLES::

                sage: from sage.misc.feature_test import SharedLibrary
                sage: emptylib = SharedLibrary("emptylib", test_code="")
                sage: emptylib.is_present()
                FeatureTestResult('emptylib', True)
        """
        # There seems to be no module readily available in python which does library checks of this kind.
        # At least other projects have also just rolled their own tests, see, e.g.:
        # http://stackoverflow.com/questions/28843765/setup-py-check-if-non-python-library-dependency-exists
        # http://stackoverflow.com/questions/28949136/checking-native-dependency-library-is-installed-in-python-setup-py
        from tempfile import NamedTemporaryFile
        with NamedTemporaryFile(suffix=".pyx") as pyx:
            pyx.write(self.test_code)
            pyx.flush()
            from sage.misc.cython import cython_import
            try:
                cython_import(pyx.name)
            except:
                return FeatureTestResult(self, False, reason="Failed to compile and import test code.")
        return FeatureTestResult(self, True, reason="Test code compiled and imported.")

class Module(Feature):
    r"""
    A :class:`Feature` which describes whether a python module can be imported.

    EXAMPLES:

    Not all builds of python include the ``ssl`` module, so you could check
    whether it is available::

        sage: from sage.misc.feature_test import Module
        sage: Module("ssl").require() # not tested - output depends on the python build
    """
    def __init__(self, name, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import LibFES, Module
            sage: isinstance(LibFES(), Module) # indirect doctest
            True
        """
        Feature.__init__(self, name, spkg=spkg, url=url)

    def is_present(self):
        r"""
        Return whether the module can be imported. This is determined by
        actually importing it.

        EXAMPLES::

            sage: from sage.misc.feature_test import Module
            sage: Module("sys").is_present()
            FeatureTestResult('sys', True)
        """
        import importlib
        try:
            importlib.import_module(self.name)
        except ImportError:
            return FeatureTestResult(self, False, reason="Failed to import `{name}`.".format(name=self.name))
        return FeatureTestResult(self, True, reason="Successfully imported `{name}`.".format(name=self.name))

class OptionalModule(Module):
    r"""
    A :class:`Feature` which describes whether a module has been enabled for
    this build of Sage and is functional, i.e., whether all its
    ``dependencies`` are present.

    EXAMPLES:

    The module :module:`sage.libs.fes` relies on the library libFES::
    However, even if the library is present at runtime, the module might not
    have been built because it was not present at built time. An
    :class:`OptionalModule` checks for both, whether the module was built and
    whether it is working, i.e., whether the library it relies on is present::

        sage: from sage.misc.feature_test import LibFES
        sage: LibFES().require() # optional: fes
    """
    def __init__(self, name, dependencies=[], spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import LibFES, OptionalModule
            sage: isinstance(LibFES(), OptionalModule)
            True
        """
        Module.__init__(self, name, spkg=spkg, url=url)
        self.dependencies = dependencies

    def is_present(self):
        r"""
        Return whether the optional module is present and its dependencies are
        satisfied. This is determined by checking the dependencies and actually
        importing the module.

        EXAMPLES::

           sage: from sage.misc.feature_test import LibFES
           sage: LibFES().is_present() # optional: fes
           FeatureTestResult('sage.libs.fes', True)
        """
        for dependency in self.dependencies:
            presence = dependency.is_present()
            if not presence:
                return FeatureTestResult(self, False, reason="Dependency `{dependency}` is not satisfied: {reason}".format(dependency=dependency.name, reason=presence.reason), resolution=presence.resolution)
        return super(OptionalModule, self).is_present()

class GapPackage(Feature):
    r"""
    A feature describing the presence of a GAP package.

    EXMAPLES::

        sage: from sage.misc.feature_test import GapPackage
        sage: GapPackage("grape", spkg="gap_packages")
        Feature('GAP package grape')
    """
    def __init__(self, package, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import GapPackage
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), GapPackage)
            True
        """
        Feature.__init__(self, "GAP package {package}".format(package=package), spkg=spkg, url=url)
        self.package = package

    @cached_method
    def is_present(self):
        r"""
        Return whether the package is available in GAP.

        This does not check whether this package is functional.

        EXAMPLES::

            sage: from sage.misc.feature_test import GapPackage
            sage: GapPackage("grape", spkg="gap_packages").is_present() # optional: gap_packages
            FeatureTestResult('GAP package grape', True)
        """
        from sage.libs.gap.libgap import libgap
        command = 'TestPackageAvailability("{package}")'.format(package=self.package)
        presence = libgap.eval(command)
        if presence:
            return FeatureTestResult(self, True,
                    reason = "`{command}` evaluated to `{presence}` in GAP.".format(command=command, presence=presence))
        else:
            return FeatureTestResult(self, False,
                    reason = "`{command}` evaluated to `{presence}` in GAP.".format(command=command, presence=presence))

class CSDP(Executable):
    r"""
    A class:`sage.misc.feature_test.Feature` which checks for the ``theta`` binary
    of CSDP.

    EXAMPLES::

        sage: from sage.misc.feature_test import CSDP
        sage: CSDP().is_present() # optional: csdp
        FeatureTestResult('CSDP', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import CSDP
            sage: isinstance(CSDP(), CSDP)
            True
        """
        Executable.__init__(self, name="CSDP", spkg="csdp", executable="theta", url="http://github.org/dimpase/csdp")

    def is_functional(self):
        r"""
        Check whether ``theta`` works on a trivial example.

        EXAMPLES::

            sage: from sage.misc.feature_test import CSDP
            sage: CSDP().is_functional() # optional: csdp
            FeatureTestResult('CSDP', True)
        """
        from sage.misc.feature_test import FeatureTestResult
        from sage.misc.temporary_file import tmp_filename
        import os, subprocess
        tf_name = tmp_filename()
        with open(tf_name, 'wb') as tf:
            tf.write("2\n1\n1 1")
        devnull = open(os.devnull, 'wb')
        command = ['theta', tf_name]
        try:
            lines = subprocess.check_output(command, stderr=devnull)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                reason = "Call to `{command}` failed with exit code {e.returncode}.".format(command=" ".join(command), e=e))

        result = lines.strip().split('\n')[-1]
        import re
        match = re.match("^The Lovasz Theta Number is (.*)$", result)
        if match is None:
            return FeatureTestResult(self, False,
                reason = "Last line of the output of `{command}` did not have the expected format.".format(command=" ".join(command)))

        return FeatureTestResult(self, True)

class Plantri(Executable):
    r"""
    A class:`sage.misc.feature_test.Feature` which checks for the ``plantri``
    binary.

    EXAMPLES::

        sage: from sage.misc.feature_test import Plantri
        sage: Plantri().is_present() # optional: plantri
        FeatureTestResult('plantri', False)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Plantri
            sage: isinstance(Plantri(), Plantri)
            True
        """
        Executable.__init__(self, name="plantri", spkg="plantri", executable="plantri", url="http://users.cecs.anu.edu.au/~bdm/plantri/")

    def is_functional(self):
        r"""
        Check whether ``plantri`` works on trivial input.

        EXAMPLES::

            sage: from sage.misc.feature_test import Plantri
            sage: Plantri().is_functional() # optional: plantri
            FeatureTestResult('plantri', False)
        """
        from sage.misc.feature_test import FeatureTestResult
        import os, subprocess
        command = ["plantri", "4"]
        try:
            lines = subprocess.check_output(command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = "1 triangulation written"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                    reason = "Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)

class Lrs(Executable):
    r"""
    A :class:`sage.misc.feature_test.Feature` describing the presence of the ``lrs``
    binary which comes as a part of ``lrslib``.

    EXAMPLES::

        sage: from sage.misc.feature_test import Lrs
        sage: Lrs().is_present() # optional: lrslib
        FeatureTestResult('lrslib', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Lrs
            sage: isinstance(Lrs(), Lrs)
            True
        """
        Executable.__init__(self, "lrslib", executable="lrs", spkg="lrslib", url="http://cgm.cs.mcgill.ca/~avis/C/lrs.html")

    def is_functional(self):
        r"""
        Test whether ``lrs`` works on a trivial input.

        EXAMPLES::

            sage: from sage.misc.feature_test import Lrs
            sage: Lrs().is_functional() # optional: lrslib
            FeatureTestResult('lrslib', True)
        """
        from sage.misc.feature_test import FeatureTestResult
        from sage.misc.temporary_file import tmp_filename
        import os, subprocess
        tf_name = tmp_filename()
        with open(tf_name, 'wb') as tf:
            tf.write("V-representation\nbegin\n 1 1 rational\n 1 \nend\nvolume")
        devnull = open(os.devnull, 'wb')
        command = ['lrs', tf_name]
        try:
            lines = subprocess.check_output(command, stderr=devnull)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                reason = "Call to `{command}` failed with exit code {e.returncode}.".format(command=" ".join(command), e=e))

        expected = "Volume= 1"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                reason = "Output of `{command}` did not contain the expected result `{expected}`.".format(command=" ".join(command),expected=expected))

        return FeatureTestResult(self, True)

class Buckygen(Executable):
    r"""
    A class:`sage.misc.feature_test.Feature` which checks for the ``buckygen``
    binary.

    EXAMPLES::

        sage: from sage.misc.feature_test import Buckygen
        sage: Buckygen().is_present() # optional: buckygen
        FeatureTestResult('Buckygen', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Buckygen
            sage: isinstance(Buckygen(), Buckygen)
            True
        """
        Executable.__init__(self, name="Buckygen", spkg="buckygen", executable="buckygen", url="http://caagt.ugent.be/buckygen/")

    def is_functional(self):
        r"""
        Check whether ``buckygen`` works on trivial input.

        EXAMPLES::

            sage: from sage.misc.feature_test import Buckygen
            sage: Buckygen().is_functional() # optional: buckygen
            FeatureTestResult('Buckygen', True)
        """
        from sage.misc.feature_test import FeatureTestResult
        import subprocess
        command = ["buckygen", "-d", "22d"]
        try:
            lines = subprocess.check_output(command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason = "Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = "Number of fullerenes generated with 13 vertices: 0"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                    reason = "Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)

class Benzene(Executable):
    r"""
    A class:`sage.misc.feature_test.Feature` which checks for the ``benzene``
    binary.

    EXAMPLES::

        sage: from sage.misc.feature_test import Benzene
        sage: Benzene().is_present() # optional: benzene
        True
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Benzene
            sage: isinstance(Benzene(), Benzene)
            True
        """
        Executable.__init__(self, name="Benzene", spkg="benzene", executable="benzene", url="http://www.grinvin.org/")

    def is_functional(self):
        r"""
        Check whether ``benzene`` works on trivial input.

        EXAMPLES::

            sage: from sage.misc.feature_test import Benzene
            sage: Benzene().is_functional() # optional: benzene
            True
        """
        from sage.misc.feature_test import FeatureTestResult
        import os, subprocess
        devnull = open(os.devnull, 'wb')
        command = ["benzene", "2", "p"]
        try:
            lines = subprocess.check_output(command, stderr=devnull)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = ">>planar_graph<<"
        if not lines.startswith(expected):
            return FeatureTestResult(self, False,
                    reason="Call `{command}` did not produce output that started with `{expected}`.".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)

class DatabaseCremona(StaticFile):
    r"""
    A :class:`Feature` which describes the presence of John Cremona's database
    of elliptic curves.

    EXAMPLES::

        sage: from sage.misc.feature_test import DatabaseCremona
        sage: DatabaseCremona().is_present() # optional: database_cremona_ellcurve
        FeatureTestResult("Cremona's database of elliptic curves", True)
    """
    def __init__(self, name="cremona", spkg=None):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import DatabaseCremona
            sage: isinstance(DatabaseCremona(), DatabaseCremona)
            True
        """
        if spkg is None and name == "cremona":
            spkg = "database_cremona_ellcurve"
        filename = name.replace(' ','_') + ".db"
        import os.path, os
        from sage.env import SAGE_SHARE
        # Ideally we would also search for the cremona db in other places but SAGE_SHARE.
        # However, the calling  code currently only looks for it there.
        StaticFile.__init__(self, "Cremona's database of elliptic curves", filename=filename, search_path=[os.path.join(SAGE_SHARE, 'cremona')], spkg=spkg, url="https://github.com/JohnCremona/ecdata")

class SmallGroupsLibrary(Feature):
    r"""
    A feature describing the presence of the Small Groups Library for GAP.

    EXMAPLES::

        sage: from sage.misc.feature_test import SmallGroupsLibrary
        sage: SmallGroupsLibrary()
        Feature('Small Groups Library')
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import SmallGroupsLibrary
            sage: isinstance(SmallGroupsLibrary(), SmallGroupsLibrary)
            True
        """
        Feature.__init__(self, "Small Groups Library", spkg="database_gap", url="www.gap-system.org/Packages/sgl.html")

    @cached_method
    def is_present(self):
        r"""
        Return whether the Small Groups Library is available in GAP.

        EXAMPLES::

            sage: from sage.misc.feature_test import SmallGroupsLibrary
            sage: SmallGroupsLibrary().is_present() # optional: database_gap
            FeatureTestResult('Small Groups Library', True)
        """
        from sage.libs.gap.libgap import libgap
        command = 'SmallGroup(13,1)'
        output = None
        presence = False
        try:
            output = str(libgap.eval(command))
            presence = True
        except ValueError as e:
            output = str(e)
        return FeatureTestResult(self, presence,
            reason = "`{command}` evaluated to `{output}` in GAP.".format(command=command, output=output))

class LibFESLibrary(SharedLibrary):
    r"""
    A :class:`Feature` which describes whether the FES library
    is present and functional.

    EXAMPLES::

        sage: from sage.misc.feature_test import LibFESLibrary
        sage: LibFESLibrary().require() # optional: fes
    """
    _test_code=r"""
#clib fes
from libc.stdint cimport uint64_t
cdef extern from "<fes_interface.h>":
    ctypedef int (*solution_callback_t)(void *, uint64_t)
    void exhaustive_search_wrapper(int n, int n_eqs, int degree, int ***coeffs, solution_callback_t callback, void* callback_state, int verbose)

solutions = 0

class InternalState:
    verbose = False
    sols = []
    max_sols = 0

cdef int report_solution(void *_state, uint64_t i):
    global solutions
    solutions += 1
    return 0

sig_on()
cdef int ***coeffs = <int ***> sig_calloc(1, sizeof(int **))
coeffs[0] = <int **> sig_calloc(3, sizeof(int *))
coeffs[0][0] = <int *> sig_calloc(1, sizeof(int))
coeffs[0][1] = <int *> sig_calloc(2, sizeof(int))
coeffs[0][2] = <int *> sig_calloc(1, sizeof(int))
coeffs[0][2][0] = 1 # x*y = 0
internal_state = InternalState()

exhaustive_search_wrapper(2, 1, 2, coeffs, report_solution, <void *>internal_state, 0)

sig_free(coeffs[0][2])
sig_free(coeffs[0][1])
sig_free(coeffs[0][0])
sig_free(coeffs[0])
sig_free(coeffs)
sig_off()

if solutions != 3: raise ImportError("libFES did not find three solutions for x*y = 0.")
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import LibFESLibrary
            sage: isinstance(LibFESLibrary(), LibFESLibrary)
            True
        """
        SharedLibrary.__init__("LibFES", test_code=LibFESLibrary._test_code, spkg="fes", url="http://www.lifl.fr/~bouillag/fes/")

class LibFES(OptionalModule):
    r"""
    A :class:`Feature` which describes whether the ;module:`sage.libs.fes`
    module has been enabled for this build of Sage and is functional.

    .. NOTE::

        This module does not depend on :class:`LibFESLibrary`. If the module
        imports we just assume that the library is working.

    EXAMPLES::

        sage: from sage.misc.feature_test import LibFES
        sage: LibFES().require() # optional: fes
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import LibFES
            sage: isinstance(LibFES(), LibFES)
            True
        """
        OptionalModule.__init__(self, "sage.libs.fes", spkg="fes", url="http://www.lifl.fr/~bouillag/fes/")

class BlissLibrary(SharedLibrary):
    r"""
    A :class:`Feature` which describes whether the Bliss library is present and
    functional.


    EXAMPLES::

        sage: from sage.misc.feature_test import BlissLibrary
        sage: BlissLibrary().require() # optional: bliss

    """
    _test_code=r"""
#clang c++
#clib bliss
cdef extern from "bliss/graph.hh" namespace "bliss":
    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)

sig_on()
Graph(1)
sig_off()
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import BlissLibrary
            sage: BlissLibrary()
            Feature('Bliss')

        """
        SharedLibrary.__init__("Bliss", test_code=Bliss._test_code, spkg="bliss", url="http://www.tcs.hut.fi/Software/bliss/")

class Bliss(OptionalModule):
    r"""
    A :class:`Feature` which describes whether the ;module:`sage.graphs.bliss`
    module has been enabled for this build of Sage and is functional.

    .. NOTE::

        This module does not depend on :class:`LibFESLibrary`. If the module
        imports we just assume that the library is working.

    EXAMPLES::

        sage: from sage.misc.feature_test import Bliss
        sage: Bliss().require() # optional: bliss

    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature_test import Bliss
            sage: Bliss()
            Feature('sage.graphs.bliss')

        """
        OptionalModule.__init__(self, "sage.graphs.bliss", spkg="bliss", url="http://www.tcs.hut.fi/Software/bliss/")
