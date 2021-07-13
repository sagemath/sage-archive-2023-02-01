r"""
Base class for mutable polyhedra.

Just like vectors and matrices they can be set immutable.
The constructor does this by default.
"""

from .base import Polyhedron_base
from sage.misc.lazy_attribute import lazy_attribute


class Polyhedron_mutable(Polyhedron_base):
    """
    Base class for polyhedra that allow mutability.

    This should not be used directly.
    """

    def __hash__(self):
        r"""
        TESTS::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: set([p])
            Traceback (most recent call last):
            ...
            TypeError: mutable polyhedra are unhashable
            sage: p.set_immutable()
            sage: set([p])
            {A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex}
        """
        if self._is_mutable:
            raise TypeError("mutable polyhedra are unhashable")
        return Polyhedron_base.__hash__(self)

    def _clear_cache(self):
        r"""
        Clear the Vrepresentation and Hrepresentation data of ``self``.

        TESTS::

            sage: p = polytopes.permutahedron(4)
            sage: P = p.parent()
            sage: q = P._element_constructor_(p, mutable=True)
            sage: TestSuite(q).run()
            sage: q._clear_cache()
            sage: TestSuite(q).run()

        ::

            sage: q.set_immutable()
            sage: q._clear_cache()
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra
        """
        if not self._is_mutable:
            raise TypeError("cannot clear cache of immutable polyhedra")
        self.parent().recycle(self)
        backend_object = self.__dict__["_" + self._backend_object_name]
        del self.__dict__
        self.__dict__["_" + self._backend_object_name] = backend_object
        self._is_mutable = True

    def is_mutable(self):
        r"""
        Return True if the polyhedron is mutable, i.e. it can be modified in place.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: p.is_mutable()
            True
            sage: p = Polyhedron([[1, 1]], mutable=False)
            sage: p.is_mutable()
            False
        """
        return self._is_mutable

    def is_immutable(self):
        r"""
        Return True if the polyhedron is immutable, i.e. it cannot be modified in place.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: p.is_immutable()
            False
            sage: p = Polyhedron([[1, 1]], mutable=False)
            sage: p.is_immutable()
            True
        """
        return not self._is_mutable

    def set_immutable(self):
        r"""
        Make this polyhedron immutable. This operation cannot be undone.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.set_immutable(p)
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        raise NotImplementedError("a derived class must implement this")

    def Vrepresentation(self):
        r"""
        A derived class must overwrite such that it restores Vrepresentation
        after clearing it.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.Vrepresentation(p)
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

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.Hrepresentation(p)
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        # A derived class must implemented it to recalculate, if necessary.
        raise NotImplementedError("a derived class must implement this")
