# -*- coding: utf-8 -*-
r"""
Check for GAP features
"""

from . import Feature, FeatureTestResult


class GapPackage(Feature):
    r"""
    A feature describing the presence of a GAP package.

    EXMAPLES::

        sage: from sage.features.gap import GapPackage
        sage: GapPackage("grape", spkg="gap_packages")
        Feature('GAP package grape')
    """
    def __init__(self, package, **kwds):
        r"""
        TESTS::

            sage: from sage.features.gap import GapPackage
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), GapPackage)
            True
        """
        Feature.__init__(self, "GAP package {package}".format(package=package), **kwds)
        self.package = package

    def _is_present(self):
        r"""
        Return whether the package is available in GAP.

        This does not check whether this package is functional.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape", spkg="gap_packages").is_present()  # optional: gap_packages
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


class SmallGroupsLibrary(Feature):
    r"""
    A feature describing the presence of the Small Groups Library for GAP.

    EXMAPLES::

        sage: from sage.features.gap import SmallGroupsLibrary
        sage: SmallGroupsLibrary()
        Feature('Small Groups Library')
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.gap import SmallGroupsLibrary
            sage: isinstance(SmallGroupsLibrary(), SmallGroupsLibrary)
            True
        """
        Feature.__init__(self, "Small Groups Library", spkg="database_gap", url="www.gap-system.org/Packages/sgl.html")

    def _is_present(self):
        r"""
        Return whether the Small Groups Library is available in GAP.

        EXAMPLES::

            sage: from sage.features.gap import SmallGroupsLibrary
            sage: SmallGroupsLibrary().is_present()  # optional: database_gap
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


class PrimitiveGroupsLibrary(Feature):
    r"""
    A feature describing the presence of the Primitive Groups Library for GAP.

    EXMAPLES::

        sage: from sage.features.gap import PrimitiveGroupsLibrary
        sage: PrimitiveGroupsLibrary()
        Feature('Primitive Groups Library')
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.gap import PrimitiveGroupsLibrary
            sage: isinstance(PrimitiveGroupsLibrary(), PrimitiveGroupsLibrary)
            True
        """
        Feature.__init__(self, "Primitive Groups Library", spkg="database_gap", url="www.gap-system.org/Datalib/prim.html")

    def _is_present(self):
        r"""
        Return whether the Primitive Groups Library is available in GAP.

        EXAMPLES::

            sage: from sage.features.gap import PrimitiveGroupsLibrary
            sage: PrimitiveGroupsLibrary().is_present()  # optional: database_gap
            FeatureTestResult('Primitive Groups Library', True)
        """
        from sage.libs.gap.libgap import libgap
        command = 'PrimitiveGroup(5,1)'
        output = None
        presence = False
        try:
            output = str(libgap.eval(command))
            presence = True
        except ValueError as e:
            output = str(e)
        return FeatureTestResult(self, presence,
            reason = "`{command}` evaluated to `{output}` in GAP.".format(command=command, output=output))
