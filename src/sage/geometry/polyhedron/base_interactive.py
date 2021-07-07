r"""
Base class for interactive and lazy polyhedra.
"""

from .base import Polyhedron_base
from sage.misc.lazy_attribute import lazy_attribute


class Polyhedron_interactive(Polyhedron_base):

    def _clear_cache(self):
        r"""
        Clear the Vrepresentation and Hrepresentation data of ``self``.

        TESTS::

            sage: P = polytopes.permutahedron(4)
            sage: TestSuite(P).run()
            sage: P._clear_cache()
            sage: TestSuite(P).run()
        """
        self.parent().recycle(self)
        backend_object = self.__dict__["_" + self._backend_object_name]
        del self.__dict__
        self.__dict__["_" + self._backend_object_name] = backend_object

    def Vrepresentation(self):
        r"""
        A derived class must overwrite such that it restores Vrepresentation
        after clearing it.

        TESTS::

            sage: from sage.geometry.polyhedron.base_interactive import Polyhedron_interactive
            sage: p = polytopes.cube()
            sage: Polyhedron_interactive.Vrepresentation(p)
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        # A derived class must implemented it to recalculate, if necessary.
        raise NotImplementedError("a derived class must implement this")

    def Hrepresentation(self):
        r"""
        A derived class must overwrite such that it restores Hrepresentation
        after clearing it.

        TESTS::

            sage: from sage.geometry.polyhedron.base_interactive import Polyhedron_interactive
            sage: p = polytopes.cube()
            sage: Polyhedron_interactive.Hrepresentation(p)
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        # A derived class must implemented it to recalculate, if necessary.
        raise NotImplementedError("a derived class must implement this")
