r"""
Base class for parent objects

CLASS HIEARCHY:

SageObject
    CategoryObject
        Parent


TESTS:
This came up in some subtle bug once.
    sage: gp(2) + gap(3)
    5
"""
cimport element
cimport sage.categories.morphism as morphism
cimport sage.categories.map as map

from copy import copy

cdef int bad_parent_warnings = 0
cdef int unique_parent_warnings = 0

# TODO: define this once?

cdef object elt_parent = None

cdef inline parent_c(x):
    if PY_TYPE_CHECK(x, element.Element):
        return (<element.Element>x)._parent
    elif hasattr(x, 'parent'):
        return x.parent()
    else:
        return <object>PY_TYPE(x)

cdef _record_exception():
    from element import get_coercion_model
    get_coercion_model()._record_exception()

cdef object _Integer
cdef bint is_Integer(x):
    global _Integer
    if _Integer is None:
        from sage.rings.integer import Integer as _Integer
    return PY_TYPE_CHECK_EXACT(x, _Integer) or PY_TYPE_CHECK_EXACT(x, int)

# for override testing
cdef extern from "descrobject.h":
    ctypedef struct PyMethodDef:
        void *ml_meth
    ctypedef struct PyMethodDescrObject:
        PyMethodDef *d_method
    void* PyCFunction_GET_FUNCTION(object)
    bint PyCFunction_Check(object)

cdef extern from *:
    Py_ssize_t PyDict_Size(object)
    Py_ssize_t PyTuple_GET_SIZE(object)

    ctypedef class __builtin__.dict [object PyDictObject]:
        cdef Py_ssize_t ma_fill
        cdef Py_ssize_t ma_used

    void* PyDict_GetItem(object, object)

cdef inline Py_ssize_t PyDict_GET_SIZE(o):
    return (<dict>o).ma_used

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

import operator
import weakref

from category_object import CategoryObject
from generators import Generators

# TODO: Theses should probably go in the printer module (but lots of stuff imports them from parent)
from category_object import localvars

cdef object BuiltinMethodType = type(repr)

include '../ext/python_object.pxi'
include '../ext/python_bool.pxi'
include '../ext/stdsage.pxi'

def is_Parent(x):
    """
    Return True if x is a parent object, i.e., derives from
    sage.structure.parent.Parent and False otherwise.

    EXAMPLES:
        sage: from sage.structure.parent import is_Parent
        sage: is_Parent(2/3)
        False
        sage: is_Parent(ZZ)
        True
        sage: is_Parent(Primes())
        True
    """
    return PY_TYPE_CHECK(x, Parent)

cdef bint guess_pass_parent(parent, element_constructor):
    if PY_TYPE_CHECK(element_constructor, types.MethodType):
        return False
    elif PY_TYPE_CHECK(element_constructor, BuiltinMethodType):
        return element_constructor.__self__ is not parent
    else:
        return True



cdef class Parent(category_object.CategoryObject):
    """
    Parents are the SAGE/mathematical analogues of container objects
    in computer science.
    """

    def __init__(self, base=None, *, categories=[], element_constructor=None, gens=None, names=None, normalize=True, **kwds):
        CategoryObject.__init__(self, categories, base)
        if len(kwds) > 0:
            if bad_parent_warnings:
                print "Illegal keywords for %s: %s" % (type(self), kwds)
        # TODO: many classes don't call this at all, but __new__ crashes SAGE
        if bad_parent_warnings:
            if element_constructor is None:
                element_constructor = self._element_constructor_
            elif not callable(element_constructor):
                print "coerce BUG: Bad element_constructor provided", type(self), type(element_constructor), element_constructor
        if gens is not None:
            self._populate_generators_(gens, names, normalize)
        elif names is not None:
            self._assign_names(names, normalize)
        self._element_constructor = element_constructor
        self._element_init_pass_parent = guess_pass_parent(self, element_constructor)
        self._coerce_from_list = []
        self._coerce_from_hash = {}
        self._action_list = []
        self._action_hash = {}
        self._convert_from_list = []
        self._convert_from_hash = {}
        self._embedding = None
        self._initial_coerce_list = []
        self._initial_action_list = []
        self._initial_convert_list = []

    cdef int init_coerce(self, bint warn=True) except -1:
        if self._coerce_from_hash is None:
            if warn:
                print "init_coerce() for ", type(self)
                raise ZeroDivisionError, "hello"
            self._initial_coerce_list = []
            self._initial_action_list = []
            self._initial_convert_list = []
            self._coerce_from_list = []
            self._coerce_from_hash = {}
            self._action_list = []
            self._action_hash = {}
            self._convert_from_list = []
            self._convert_from_hash = {}
            self._embedding = None

    def _introspect_coerce(self):
        """
        Used for debugging the coercion model.

        EXAMPLES:
            sage: sorted(QQ._introspect_coerce().items())
            [('_action_hash', {...}),
             ('_action_list', []),
             ('_coerce_from_hash', {...}),
             ('_coerce_from_list', []),
             ('_convert_from_hash', {...}),
             ('_convert_from_list', [...]),
             ('_element_init_pass_parent', False),
             ('_embedding', None),
             ('_initial_action_list', []),
             ('_initial_coerce_list', []),
             ('_initial_convert_list', [])]
        """
        return {
            '_coerce_from_list': self._coerce_from_list,
            '_coerce_from_hash': self._coerce_from_hash,
            '_action_list': self._action_list,
            '_action_hash': self._action_hash,
            '_convert_from_list': self._convert_from_list,
            '_convert_from_hash': self._convert_from_hash,
            '_embedding': self._embedding,
            '_initial_coerce_list': self._initial_coerce_list,
            '_initial_action_list': self._initial_action_list,
            '_initial_convert_list': self._initial_convert_list,
            '_element_init_pass_parent': self._element_init_pass_parent,
        }

    def __getstate__(self):
        """
        Used for pickling.

        TESTS:
            sage: loads(dumps(RR['x'])) == RR['x']
            True
        """
        d = CategoryObject.__getstate__(self)
        d['_embedding'] = self._embedding
        d['_element_constructor'] = self._element_constructor
        d['_convert_method_name'] = self._convert_method_name
        d['_element_init_pass_parent'] = self._element_init_pass_parent
        d['_initial_coerce_list'] = self._initial_coerce_list
        d['_initial_action_list'] = self._initial_action_list
        d['_initial_convert_list'] = self._initial_convert_list
        return d

    def __setstate__(self, d):
        """
        Used for pickling.

        TESTS:
            sage: loads(dumps(CDF['x'])) == CDF['x']
            True
        """
        CategoryObject.__setstate__(self, d)
        try:
            version = d['_pickle_version']
        except KeyError:
            version = 0
        if version == 1:
            self.init_coerce(False) # Really, do we want to init this with the same initial data as before?
            self._populate_coercion_lists_(coerce_list=d['_initial_coerce_list'] or [],
                                           action_list=d['_initial_action_list'] or [],
                                           convert_list=d['_initial_convert_list'] or [],
                                           embedding=d['_embedding'],
                                           convert_method_name=d['_convert_method_name'],
                                           element_constructor = d['_element_constructor'],
                                           init_no_parent=not d['_element_init_pass_parent'],
                                           unpickling=True)

    def __call__(self, x=0, *args, **kwds):
        """
        This is the generic call method for all parents.

        When called, it will find a map based on the Parent (or type) of x.
        If a coercion exists, it will always be chosen. This map will
        then be called (with the arguments and keywords if any).

        By default this will dispatch as quickly as possible to
        \code{self._element_constructor_()} though faster pathways are
        possible if so desired.
        """
        if self._element_constructor is None:
            raise NotImplementedError
        cdef Py_ssize_t i
        cdef R = parent_c(x)
        cdef bint no_extra_args = PyTuple_GET_SIZE(args) == 0 and PyDict_GET_SIZE(kwds) == 0
        if R is self and no_extra_args:
            return x

        # Here we inline the first part of convert_map_from for speed.
        # (Yes, the virtual function overhead can matter.)
        if self._convert_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce()
        cdef map.Map mor
        cdef PyObject* mor_ptr = PyDict_GetItem(self._convert_from_hash, R)
        if mor_ptr != NULL:
            mor = <map.Map>mor_ptr
        else:
            mor = <map.Map>self.convert_map_from(R)

        if mor is not None:
            if no_extra_args:
                return mor._call_(x)
            else:
                return mor._call_with_args(x, args, kwds)

        raise TypeError, "No conversion defined from %s to %s"%(R, self)

    #############################################################################
    # Containment testing
    #############################################################################
    def __contains__(self, x):
        r"""
        True if there is an element of self that is equal to x under ==,
        or if x is already an element of self.

        For many structures we test this by using \code{__call__} and
        then testing equality between x and the result.

        EXAMPLES:
            sage: 2 in Integers(7)
            True
            sage: 2 in ZZ
            True
            sage: Integers(7)(3) in ZZ
            True
            sage: 3/1 in ZZ
            True
            sage: 5 in QQ
            True
            sage: I in RR
            False
            sage: SR(2) in ZZ
            True
            sage: RIF(1, 2) in RIF
            True
            sage: pi in RIF # there is no element of RIF equal to pi
            False
        """
        if parent_c(x) == self:
            return True
        try:
            x2 = self(x)
            return bool(x2 == x)
        except (TypeError, ValueError):
            return False

    cpdef coerce(self, x):
        """
        Return x as an element of self, if and only if there is a canonical
        coercion from the parent of x to self.

        EXAMPLES:
            sage: QQ.coerce(ZZ(2))
            2
            sage: ZZ.coerce(QQ(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Integer Ring

        We make an exception for zero:
            sage: V = GF(7)^7
            sage: V.coerce(0)
            (0, 0, 0, 0, 0, 0, 0)
        """
        mor = self.coerce_map_from(parent_c(x))
        if mor is None:
            if is_Integer(x) and not x:
                try:
                    return self(0)
                except:
                    _record_exception()
            raise TypeError, "no canonical coercion from %s to %s" % (parent_c(x), self)
        else:
            return (<map.Map>mor)._call_(x)

    def list(self):
        """
        Return a list of all elements in this object, if possible (the
        object must define an iterator).
        """
        try:
            if self._list is not None:
                return self._list
        except AttributeError:
            pass
        the_list = list(self.__iter__())
        try:
            self._list = the_list
        except AttributeError:
            pass
        return the_list

    def __nonzero__(self):
        """
        By default, all Parents are treated as True when used in an if
        statement. Override this method if other behavior is desired
        (for example, for empty sets).

        EXAMPLES:
            sage: if ZZ: print "Yes"
            Yes
        """
        return True

    cpdef bint _richcmp_helper(left, right, int op) except -2:
        """
        Compare left and right.

        EXAMPLES:
            sage: ZZ < QQ
            True
        """
        cdef int r

        if not PY_TYPE_CHECK(right, Parent) or not PY_TYPE_CHECK(left, Parent):
            # One is not a parent -- use arbitrary ordering
            if (<PyObject*>left) < (<PyObject*>right):
                r = -1
            elif (<PyObject*>left) > (<PyObject*>right):
                r = 1
            else:
                r = 0

        else:
            # Both are parents -- but need *not* have the same type.
            if HAS_DICTIONARY(left):
                r = left.__cmp__(right)
            else:
                r = left._cmp_(right)

        if op == 0:  #<
            return r  < 0
        elif op == 2: #==
            return r == 0
        elif op == 4: #>
            return r  > 0
        elif op == 1: #<=
            return r <= 0
        elif op == 3: #!=
            return r != 0
        elif op == 5: #>=
            return r >= 0

    def __len__(self):
        """
        Returns the number of elements in self. This is the naive algorithm
        of listing self and counting the elements.

        EXAMPLES:
            sage: len(GF(5))
            5
            sage: len(MatrixSpace(GF(2), 3, 3))
            512
        """
        return len(self.list())

    def __getitem__(self, n):
        """
        Returns the $n$-th item of self, by getting self as a list.

        EXAMPLES:
            sage: MatrixSpace(GF(3), 2, 2)[9]
            [0 2]
            [0 0]
            sage: MatrixSpace(GF(3), 2, 2)[0]
            [0 0]
            [0 0]
        """
        return self.list()[n]

    def __getslice__(self,  Py_ssize_t n,  Py_ssize_t m):
        """
        Returns the $n$-th through $m$-th items in self.

        EXAMPLES:
            sage: VectorSpace(GF(7), 3)[:10]
            [(0, 0, 0),
             (1, 0, 0),
             (2, 0, 0),
             (3, 0, 0),
             (4, 0, 0),
             (5, 0, 0),
             (6, 0, 0),
             (0, 1, 0),
             (1, 1, 0),
             (2, 1, 0)]
        """
        return self.list()[int(n):int(m)]

    #################################################################################
    # Generators and Homomorphisms
    #################################################################################

    def _is_valid_homomorphism_(self, codomain, im_gens):
       r"""
       Return True if \code{im_gens} defines a valid homomorphism
       from self to codomain; otherwise return False.

       If determining whether or not a homomorphism is valid has not
       been implemented for this ring, then a NotImplementedError exception
       is raised.
       """
       raise NotImplementedError, "Verification of correctness of homomorphisms from %s not yet implmented."%self

    def Hom(self, codomain, cat=None):
        r"""
        self.Hom(codomain, cat=None):

        Return the homspace \code{Hom(self, codomain, cat)} of all
        homomorphisms from self to codomain in the category cat.  The
        default category is \code{self.category()}.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: R.Hom(QQ)
            Set of Homomorphisms from Multivariate Polynomial Ring in x, y over Rational Field to Rational Field

        Homspaces are defined for very general \sage objects, even elements of familiar rings.
            sage: n = 5; Hom(n,7)
            Set of Morphisms from 5 to 7 in Category of elements of Integer Ring
            sage: z=(2/3); Hom(z,8/1)
            Set of Morphisms from 2/3 to 8 in Category of elements of Rational Field

        This example illustrates the optional third argument:
            sage: QQ.Hom(ZZ, Sets())
            Set of Morphisms from Rational Field to Integer Ring in Category of sets
        """
        try:
            return self._Hom_(codomain, cat)
        except (TypeError, AttributeError):
            pass
        from sage.categories.all import Hom
        return Hom(self, codomain, cat)

    def hom(self, im_gens, codomain=None, check=None):
       r"""
       Return the unique homomorphism from self to codomain that
       sends \code{self.gens()} to the entries of \code{im_gens}.
       Raises a TypeError if there is no such homomorphism.

       INPUT:
           im_gens -- the images in the codomain of the generators of
                      this object under the homomorphism
           codomain -- the codomain of the homomorphism
           check -- whether to verify that the images of generators extend
                    to define a map (using only canonical coercisions).

       OUTPUT:
           a homomorphism self --> codomain

       \note{As a shortcut, one can also give an object X instead of
       \code{im_gens}, in which case return the (if it exists)
       natural map to X.}

       EXAMPLE: Polynomial Ring
       We first illustrate construction of a few homomorphisms
       involving a polynomial ring.

           sage: R.<x> = PolynomialRing(ZZ)
           sage: f = R.hom([5], QQ)
           sage: f(x^2 - 19)
           6

           sage: R.<x> = PolynomialRing(QQ)
           sage: f = R.hom([5], GF(7))
           Traceback (most recent call last):
           ...
           TypeError: images do not define a valid homomorphism

           sage: R.<x> = PolynomialRing(GF(7))
           sage: f = R.hom([3], GF(49,'a'))
           sage: f
           Ring morphism:
             From: Univariate Polynomial Ring in x over Finite Field of size 7
             To:   Finite Field in a of size 7^2
             Defn: x |--> 3
           sage: f(x+6)
           2
           sage: f(x^2+1)
           3

       EXAMPLE: Natural morphism
           sage: f = ZZ.hom(GF(5))
           sage: f(7)
           2
           sage: f
           Ring Coercion morphism:
             From: Integer Ring
             To:   Finite Field of size 5

       There might not be a natural morphism, in which case a TypeError exception is raised.
           sage: QQ.hom(ZZ)
           Traceback (most recent call last):
           ...
           TypeError: Natural coercion morphism from Rational Field to Integer Ring not defined.
       """
       if isinstance(im_gens, Parent):
           return self.Hom(im_gens).natural_map()
       from sage.structure.all import Sequence
       if codomain is None:
           im_gens = Sequence(im_gens)
           codomain = im_gens.universe()
       if isinstance(im_gens, (Sequence, Generators)):
            im_gens = list(im_gens)
       if check is None:
           return self.Hom(codomain)(im_gens)
       else:
           return self.Hom(codomain)(im_gens, check=check)

    #################################################################################
    # New New Coercion support functionality
    #################################################################################

    def _populate_coercion_lists_(self,
                                  coerce_list=[],
                                  action_list=[],
                                  convert_list=[],
                                  embedding=None,
                                  convert_method_name=None,
                                  element_constructor=None,
                                  init_no_parent=None,
                                  bint unpickling=False):
        """
        This function allows one to specify coercions, actions, conversions
        and embeddings involving this parent.

        IT SHOULD ONLY BE CALLED DURING THE __INIT__ method, often at the end.

        INPUT:
            coerce_list    -- a list of coercion Morphisms to self and
                              parents with canonical coercions to self
            action_list    -- a list of actions on and by self
            convert_list   -- a list of conversion Maps to self and
                              parents with conversions to self
            embedding      -- a single Morphism from self
            convert_method_name -- a name to look for that other elements
                              can implement to create elements of self (e.g. _integer_)
            element_constructor -- A callable object used by the __call__ method to
                              construct new elements. Typically the element class or a
                              bound method (defaults to self._element_constructor_).
            init_no_parent -- if True omit passing self in as the first
                              argument of element_constructor for conversion. This is
                              useful if parents are unique, or element_constructor is
                              a bound method (this latter case can be detected
                              automatically).
        """
        self.init_coerce(False)

        if element_constructor is None and not unpickling:
            try:
                element_constructor = self._element_constructor_
            except AttributeError:
                raise RuntimeError, "element_constructor must be provided, either as an _element_constructor_ method or via the _populate_coercion_lists_ call"
        self._element_constructor = element_constructor
        self._element_init_pass_parent = guess_pass_parent(self, element_constructor)

        if not PY_TYPE_CHECK(coerce_list, list):
            raise ValueError, "%s_populate_coercion_lists_: coerce_list is type %s, must be list" % (type(coerce_list), type(self))
        if not PY_TYPE_CHECK(action_list, list):
            raise ValueError, "%s_populate_coercion_lists_: action_list is type %s, must be list" % (type(action_list), type(self))
        if not PY_TYPE_CHECK(convert_list, list):
            raise ValueError, "%s_populate_coercion_lists_: convert_list is type %s, must be list" % (type(convert_list), type(self))

        self._initial_coerce_list = copy(coerce_list)
        self._initial_action_list = copy(action_list)
        self._initial_convert_list = copy(convert_list)

        self._convert_method_name = convert_method_name
        if init_no_parent is not None:
            self._element_init_pass_parent = not init_no_parent

        for mor in coerce_list:
            if PY_TYPE_CHECK(mor, map.Map):
                if mor.codomain() is not self:
                    raise ValueError, "Map's codomain must be self (%s) is not (%s)" % (self, mor.codomain())
                self._coerce_from_list.append(mor)
                self._coerce_from_hash[mor.domain()] = mor
            elif PY_TYPE_CHECK(mor, Parent) or PY_TYPE_CHECK(mor, type):
                P = mor
                mor = self._generic_convert_map(mor)
                self._coerce_from_list.append(mor)
                self._coerce_from_hash[P] = mor
            else:
                raise TypeError, "entries in the coerce_list must be parents or maps (got %s)" % type(mor)

        from sage.categories.action import Action
        for action in action_list:
            if isinstance(action, Action):
                if action.actor() is self:
                    self._action_list.append(action)
                    self._action_list[action.domain(), action.operation(), action.is_left()] = action
                elif action.domain() is self:
                    self._action_list.append(action)
                    self._action_list[action.actor(), action.operation(), not action.is_left()] = action
                else:
                    raise ValueError, "Action must involve self"
            else:
                raise TypeError, "entries in the action_list must be actions"

        for mor in convert_list:
            if isinstance(mor, map.Map):
                if mor.codomain() is not self:
                    raise ValueError, "Map's codomain must be self"
                self._convert_from_list.append(mor)
                self._convert_from_hash[mor.domain()] = mor
            elif PY_TYPE_CHECK(mor, Parent) or PY_TYPE_CHECK(mor, type):
                t = mor
                mor = self._generic_convert_map(mor)
                self._convert_from_list.append(mor)
                self._convert_from_hash[t] = mor
                self._convert_from_hash[mor.domain()] = mor
            else:
                raise TypeError, "entries in the convert_list must be parents or maps"

        if isinstance(embedding, map.Map):
            if embedding.domain() is not self:
                raise ValueError, "embedding's domain must be self"
            self._embedding = embedding
        elif isinstance(embedding, Parent):
            self._embedding = embedding._generic_convert_map(self)
        elif embedding is not None:
            raise TypeError, "embedding must be a parent or map"

    def coerce_embedding(self):
        """
        Returns the embedding of self into some other parent, if such a parent
        exists.

        This does not mean that there are no coercion maps from self into other
        fields, this is simply a specific morphism specified out of self and
        ususally denotes a special relationship (e.g. sub-objects, choice of
        completion, etc.)

        EXAMPLES:
        """
        return self._embedding

    cpdef _generic_convert_map(self, S):
        r"""
        Returns the default conversion map based on the data provided to
        \code{_populate_coercion_lists_()}.

        This called when \code{_coerce_map_from_()} returns \code{True}.

        If a \code{convert_method_name} is provided, it creates a \code{NamedConvertMap}, otherwise
        it creates a \code{DefaultConvertMap} or \code{DefaultConvertMap_unique} depending
        on whether or not init_no_parent is set.

        EXAMPLES:
        """
        import coerce_maps
        if self._convert_method_name is not None:
            # handle methods like _integer_
            if PY_TYPE_CHECK(S, type):
                element_constructor = S
            elif PY_TYPE_CHECK(S, Parent):
                element_constructor = (<Parent>S)._element_constructor
                if not PY_TYPE_CHECK(element_constructor, type):
                    # if element_constructor is not an actual class, get the element class
                    element_constructor = type(S.an_element())
            else:
                element_constructor = None
            if element_constructor is not None and hasattr(element_constructor, self._convert_method_name):
                return coerce_maps.NamedConvertMap(S, self, self._convert_method_name)

        if self._element_init_pass_parent:
            return coerce_maps.DefaultConvertMap(S, self)
        else:
            return coerce_maps.DefaultConvertMap_unique(S, self)

    def _coerce_map_via(self, v, S):
        """
        This attempts to construct a morphism from S to self by passing through
        one of the items in v (tried in order).

        S may appear in the list, in which case algorithm will never progress
        beyond that point.

        This is similar in spirit to the old {{{_coerce_try}}}, and useful when
        defining _coerce_map_from_

        INPUT:
            v -- A list (iterator) of parents with coercions into self. There
                 MUST be maps provided from each item in the list to self.
            S -- the starting parent

        EXAMPLES:
            sage: CDF._coerce_map_via([ZZ, RR, CC], int)
            Composite map:
              From: Set of Python objects of type 'int'
              To:   Complex Double Field
              Defn:   Native morphism:
                      From: Set of Python objects of type 'int'
                      To:   Integer Ring
                    then
                      Native morphism:
                      From: Integer Ring
                      To:   Complex Double Field

            sage: CDF._coerce_map_via([ZZ, RR, CC], QQ)
            Composite map:
              From: Rational Field
              To:   Complex Double Field
              Defn:   Generic map:
                      From: Rational Field
                      To:   Real Field with 53 bits of precision
                    then
                      Native morphism:
                      From: Real Field with 53 bits of precision
                      To:   Complex Double Field

            sage: CDF._coerce_map_via([ZZ, RR, CC], CC)
            Generic map:
              From: Complex Field with 53 bits of precision
              To:   Complex Double Field
        """
        cdef Parent R
        for R in v:
            if R is S:
                return self.coerce_map_from(R)
            connecting = R.coerce_map_from(S)
            if connecting is not None:
                return self.coerce_map_from(R) * connecting

    cpdef bint has_coerce_map_from(self, S) except -2:
        """
        Return True if there is a natural map from S to self.
        Otherwise, return False.

        EXAMPLES:
            sage: RDF.has_coerce_map_from(QQ)
            True
            sage: RDF.has_coerce_map_from(QQ['x'])
            False
            sage: RDF['x'].has_coerce_map_from(QQ['x'])
            True
            sage: RDF['x,y'].has_coerce_map_from(QQ['x'])
            True
        """
        if S is self:
            return True
        elif S == self:
            if unique_parent_warnings:
                print "Warning: non-unique parents %s"%(type(S))
            return True
        return self.coerce_map_from(S) is not None

    cpdef _coerce_map_from_(self, S):
        """
        Override this method to specify coercions beyond those specified
        in coerce_list.

        If no such coercion exists, return None or False. Otherwise, it may
        return either an actual Map to use for the coercion, a callable
        (in which case it will be wrapped in a Map), or True (in which case
        a generic map will be provided).
        """
        return None

    cpdef coerce_map_from(self, S):
        """
        This returns a Map object to coerce from S to self if one exists,
        or None if no such coercion exists.

        EXAMPLES:
            sage: ZZ.coerce_map_from(int)
            Native morphism:
              From: Set of Python objects of type 'int'
              To:   Integer Ring
            sage: QQ.coerce_map_from(ZZ)
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        cdef map.Map mor
        if S is self:
            from sage.categories.homset import Hom
            return Hom(self, self).identity()

        elif isinstance(S, Set_PythonType_class):
            return self.coerce_map_from(S._type)
        if self._coerce_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce()
        cdef object ret
        try:
            return self._coerce_from_hash[S]
        except KeyError:
            pass

        if S == self:
            # non-unique parents
            if unique_parent_warnings:
                print "Warning: non-unique parents %s"%(type(S))
            return self._generic_convert_map(S)

        try:
            _register_pair(self, S)
            mor = self.discover_coerce_map_from(S)
            #if mor is not None:
            #    # Need to check that this morphism doesn't connect previously unconnected parts of the coercion diagram
            #    if self._embedding is not None and not self._embedding.codomain().has_coerce_map_from(S):
            #        # The following if statement may call this function with self and S.  If so, we want to return None,
            #        # so that it doesn't use this path for the existence of a coercion path.
            #        # We disable this for now because it is too strict
            #        pass
            #        # print "embed problem: the following morphisms connect unconnected portions of the coercion graph\n%s\n%s"%(self._embedding, mor)
            #        # mor = None
            if mor is not None:
                # NOTE: this line is what makes the coercion detection stateful
                # self._coerce_from_list.append(mor)
                pass
            self._coerce_from_hash[S] = mor
            _unregister_pair(self, S)
            return mor
        except NotImplementedError:
            _record_exception()
            return None

    cdef discover_coerce_map_from(self, S):
        """
        Precedence for discovering a coercion S -> self goes as follows:

            1. If S has an embedding into T, look for T -> self and return composition
            2. If self._coerce_map_from_(S) is NOT exactly one of
                    - DefaultConvertMap
                    - DefaultConvertMap_unique
                    - NamedConvertMap
                return this map
            3. Traverse the coercion lists looking for another map
               returning the map from step (2) if none is found.

        In the future, multiple paths may be discovered and compared.
        """
        best_mor = None
        if PY_TYPE_CHECK(S, Parent) and (<Parent>S)._embedding is not None:
            if (<Parent>S)._embedding.codomain() is self:
                return (<Parent>S)._embedding
            connecting = self.coerce_map_from((<Parent>S)._embedding.codomain())
            if connecting is not None:
                return (<Parent>S)._embedding.post_compose(connecting)

        cdef map.Map mor
        user_provided_mor = self._coerce_map_from_(S)

        if user_provided_mor is False:
            user_provided_mor = None

        elif user_provided_mor is not None:

            from sage.categories.map import Map
            from coerce_maps import DefaultConvertMap, DefaultConvertMap_unique, NamedConvertMap, CallableConvertMap

            if user_provided_mor is True:
                mor = self._generic_convert_map(S)
            elif PY_TYPE_CHECK(user_provided_mor, Map):
                mor = <map.Map>user_provided_mor
            elif callable(user_provided_mor):
                mor = CallableConvertMap(user_provided_mor)
            else:
                raise TypeError, "_coerce_map_from_ must return None, a boolean, a callable, or an explicit Map (called on %s, got %s)" % (type(self), type(user_provided_mor))

            if (PY_TYPE_CHECK_EXACT(mor, DefaultConvertMap) or
                  PY_TYPE_CHECK_EXACT(mor, DefaultConvertMap_unique) or
                  PY_TYPE_CHECK_EXACT(mor, NamedConvertMap)) and not mor._force_use:
                # If there is something better in the list, try to return that instead
                # This is so, for example, _coerce_map_from_ can return True but still
                # take advantage of the _populate_coercion_lists_ data.
                best_mor = mor
            else:
                return mor

        from sage.categories.homset import Hom

        cdef int num_paths = 1 # this is the number of paths we find before settling on the best (the one with lowest coerce_cost).
                               # setting this to 1 will make it return the first path found.
        cdef int mor_found = 0
        cdef Parent R
        # Recurse.  Note that if S is the domain of one of the maps in self._coerce_from_list,
        # we will have stuck the map into _coerce_map_hash and thus returned it already.
        for mor in self._coerce_from_list:
            if mor._domain is S:
                if best_mor is None or mor._coerce_cost < best_mor._coerce_cost:
                    best_mor = mor
                mor_found += 1
                if mor_found  >= num_paths:
                    return best_mor
            connecting = mor._domain.coerce_map_from(S)
            if connecting is not None:
                mor = mor * connecting
                if best_mor is None or mor._coerce_cost < best_mor._coerce_cost:
                    best_mor = mor
                mor_found += 1
                if mor_found  >= num_paths:
                    return best_mor

        return best_mor


    cpdef convert_map_from(self, S):
        """
        This function returns a Map from S to self, which may or may not
        succeed on all inputs. If a coercion map from S to self exists,
        then the it will be returned. If a coercion from self to S exists,
        then it will attempt to return a section of that map.

        Under the new coercion model, this is the fastest way to convert
        elements of S to elements of self (short of manually constructing
        the elements) and is used by __call__.

        EXAMPLES:
            sage: m = ZZ.convert_map_from(QQ)
            sage: m(-35/7)
            -5
            sage: parent(m(-35/7))
            Integer Ring
        """
        if self._convert_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce()
        try:
            return self._convert_from_hash[S]
        except KeyError:
            mor = self.discover_convert_map_from(S)
            self._convert_from_list.append(mor)
            self._convert_from_hash[S] = mor
            return mor

    cdef discover_convert_map_from(self, S):

        cdef map.Map mor = self.coerce_map_from(S)
        if mor is not None:
            return mor

        try:
            if PY_TYPE_CHECK(S, Parent):
                mor = S.coerce_map_from(self)
                if mor is not None:
                    mor = mor.section()
                    if mor is not None:
                        return mor
        except NotImplementedError:
            pass


        user_provided_mor = self._convert_map_from_(S)

        if user_provided_mor is not None:
            if PY_TYPE_CHECK(user_provided_mor, map.Map):
                return user_provided_mor
            elif callable(user_provided_mor):
                from coerce_maps import CallableConvertMap
                return CallableConvertMap(user_provided_mor)
            else:
                raise TypeError, "_convert_map_from_ must return a map or callable (called on %s, got %s)" % (type(self), type(user_provided_mor))

        if not PY_TYPE_CHECK(S, type) and not PY_TYPE_CHECK(S, Parent):
            # Sequences is used as a category and a "Parent"
            from sage.categories.category_types import Sequences
            from sage.structure.coerce_maps import ListMorphism
            if isinstance(S, Sequences):
                return ListMorphism(S, self.convert_map_from(list))

        mor = self._generic_convert_map(S)
        return mor

    cpdef _convert_map_from_(self, S):
        """
        Override this method to provide additional conversions beyond those
        given in convert_list.

        This function is called after coercions are attempted. If there is a
        coercion morphism in the opposite direction, one should consider
        adding a section method to that.

        This MUST return a Map from S to self, or None. If None is returned
        then a generic map will be provided.
        """
        return None

    cpdef get_action(self, S, op=operator.mul, bint self_on_left=True):
        """
        Returns an action of self on S or S on self.

        To provide additional actions, override _get_action_.

        TESTS:
            sage: M = QQ['y']^3
            sage: M.get_action(ZZ['x']['y'])
            Right scalar multiplication by Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring on Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in y over Rational Field
            sage: M.get_action(ZZ['x']) # should be None
        """
        try:
            if self._action_hash is None: # this is because parent.__init__() does not always get called
                self.init_coerce()
            return self._action_hash[S, op, self_on_left]
        except KeyError:
            pass

        action = self._get_action_(S, op, self_on_left)
        if action is None:
            action = self.discover_action(S, op, self_on_left)

        if action is not None:
            from sage.categories.action import Action
            if not isinstance(action, Action):
                raise TypeError, "get_action_impl must return None or an Action"
            # We do NOT add to the list, as this would lead to errors as in
            # the example above.

        self._action_hash[S, op, self_on_left] = action
        return action


    cdef discover_action(self, S, op, bint self_on_left):
        # G acts on S, G -> G', R -> S => G' acts on R (?)
        # NO! ZZ[x,y] acts on Matrices(ZZ[x]) but ZZ[y] does not.
        # What may be true is that if the action's desination is S, then this can be allowed.
        from sage.categories.action import Action, PrecomposedAction
        from sage.categories.homset import Hom
        from coerce_actions import LeftModuleAction, RightModuleAction
        cdef Parent R
        for action in self._action_list:
            if PY_TYPE_CHECK(action, Action) and action.operation() is op:
                if self_on_left:
                    if action.left_domain() is not self: continue
                    R = action.right_domain()
                else:
                    if action.right_domain() is not self: continue
                    R = action.left_domain()
            elif op is operator.mul:
                try:
                    R = action
                    _register_pair(x,y) # to kill circular recursion
                    if self_on_left:
                        action = LeftModuleAction(R, self) # self is acted on from right
                    else:
                        action = RightModuleAction(R, self) # self is acted on from left
                    _unregister_pair(x,y)
                    ## The following two lines are diabled to prevent the following from working:
                    ## sage: x, y = var('x,y')
                    ## sage: parent(ZZ[x][y](1)*vector(QQ[y],[1,2]))
                    ## sage: parent(ZZ[x](1)*vector(QQ[y],[1,2]))
                    ## We will hopefully come up with a way to reinsert them, because they increase the scope
                    ## of discovered actions.
                    #i = self._action_list.index(R)
                    #self._action_list[i] = action
                except TypeError:
                    _record_exception()
                    _unregister_pair(x,y)
                    continue
            else:
                continue # only try mul if not specified
            if R is S:
                return action
            else:
                connecting = R.coerce_map_from(S) # S -> R
                if connecting is not None:
                    if self_on_left:
                        return PrecomposedAction(action, None, connecting)
                    else:
                        return PrecomposedAction(action, connecting, None)

        if op is operator.mul and PY_TYPE_CHECK(S, Parent):
            from coerce_actions import LeftModuleAction, RightModuleAction, LAction, RAction
            # Actors define _l_action_ and _r_action_
            # Acted-on elements define _lmul_ and _rmul_

            # TODO: if _xmul_/_x_action_ code does stuff like
            # if self == 0:
            #    return self
            # then an_element() == 0 could be very bad.
            #

            x = self.an_element()
            y = (<Parent>S).an_element()
            #print "looking for action ", self, "<--->", S
            #print "looking for action ", x, "<--->", y

            _register_pair(x,y) # this is to avoid possible infinite loops
            if self_on_left:
                try:
                    #print "RightModuleAction", S, self
                    action = RightModuleAction(S, self) # this will test _lmul_
                    _unregister_pair(x,y)
                    #print "got", action
                    return action
                except (NotImplementedError, TypeError, AttributeError, ValueError), ex:
                    _record_exception()

                try:
                    #print "LAction"
                    z = x._l_action_(y)
                    _unregister_pair(x,y)
                    #print "got", z
                    return LAction(self, S)
                except (NotImplementedError, TypeError, AttributeError, ValueError):
                    _record_exception()

            else:
                try:
                    #print "LeftModuleAction"
                    action = LeftModuleAction(S, self) # this will test _rmul_
                    _unregister_pair(x,y)
                    #print "got", action
                    return action
                except (NotImplementedError, TypeError, AttributeError, ValueError), ex:
                    _record_exception()

                try:
                    #print "RAction"
                    z = x._r_action_(y)
                    _unregister_pair(x,y)
                    #print "got", z
                    return RAction(self, S)
                except (NotImplementedError, TypeError, AttributeError, ValueError):
                    _record_exception()


            try:
                # maybe there is a more clever way of detecting ZZ than importing here...
                from sage.rings.integer_ring import ZZ
                from sage.rings.ring import Ring
                if S is ZZ and not self.has_coerce_map_from(ZZ):
                    #print "IntegerMulAction"
                    from sage.structure.coerce_actions import IntegerMulAction
                    action = IntegerMulAction(S, self, not self_on_left)
                    #print "got", action
                    _unregister_pair(x,y)
                    return action
            except (NotImplementedError, TypeError, AttributeError, ValueError):
                _record_exception()

            _unregister_pair(x,y)

            #print "found nothing"


    cpdef _get_action_(self, S, op, bint self_on_left):
        """
        Override this method to provide an action of self on S or S on self
        beyond what was specified in action_list.

        This must return an action which accepts an element of self and an
        element of S (in the order specified by self_on_left).
        """
        return None

    def construction(self):
        """
        Returns a pair (functor, parent) such that functor(parent) return self.
        If this ring does not have a functorial construction, return None.

        EXAMPLES:
            sage: QQ.construction()
            (FractionField, Integer Ring)
            sage: f, R = QQ['x'].construction()
            sage: f
            Poly[x]
            sage: R
            Rational Field
            sage: f(R)
            Univariate Polynomial Ring in x over Rational Field
        """
        return None

    cpdef an_element(self):
        r"""
        Implementation of a function that returns an element (often non-trivial)
        of a parent object.  This is cached. Parent structures that are should
        override \code{_an_element_} instead.

        EXAMPLES:
            sage: CDF.an_element()
            1.0*I
            sage: ZZ[['t']].an_element()
            t
        """
        if self.__an_element is None:
            self.__an_element = self._an_element_()
        return self.__an_element

    cpdef _an_element_(self):
        """
        Returns an element of self. Want it in sufficent generality
        that poorly-written functions won't work when they're not
        supposed to. This is cached so doesn't have to be super fast.

        EXAMPLES:
            sage: QQ._an_element_()
            1/2
            sage: ZZ['x,y,z']._an_element_()
            x
        """
        try:
            return self.gen(0)
        except:
            _record_exception()
            pass

        try:
            return self.gen()
        except:
            _record_exception()
            pass

        from sage.rings.infinity import infinity
        for x in ['_an_element_', 'pi', 1.2, 2, 1, 0, infinity]:
            # This weird looking list is to try to get an element
            # which doesn't coerce other places.
            try:
                return self(x)
            except (TypeError, NameError, NotImplementedError, AttributeError, ValueError):
                _record_exception()

        raise NotImplementedError, "please implement _an_element_ for %s" % self

    cpdef bint is_exact(self) except -2:
        """
        Return True if elements of this ring are represented exactly, i.e.,
        there is no precision loss when doing arithmetic.

        NOTE: This defaults to true, so even if it does return True you have
        no guarantee (unless the ring has properly overloaded this).

        EXAMPLES:
            sage: QQ.is_exact()
            True
            sage: ZZ.is_exact()
            True
            sage: Qp(7).is_exact()
            False
            sage: Zp(7, type='capped-abs').is_exact()
            False
        """
        return True

#    cpdef base_extend(self, other, category=None):
#        """
#        EXAMPLES:
#            sage: QQ.base_extend(GF(7))
#            Traceback (most recent call last):
#            ...
#            TypeError: base extension not defined for Rational Field
#            sage: ZZ.base_extend(GF(7))
#            Finite Field of size 7
#        """
#        # Use the coerce map if a base extend is not defined in the category.
#        # this is the default implementation.
#        try:
#            if category is None:
#                method = self._categories[0].get_object_method("base_extend") # , self._categories[1:])
#            else:
#                method = category.get_object_method("base_extend")
#            if method is not None:
#                return method(self)
#            elif other.has_coerce_map_from(self):
#                return other
#            else:
#                raise TypeError, "base extension not defined for %s" % self
#        except AttributeError:
#            raise TypeError, "base extension not defined for %s" % self

############################################################################
# Set baseclass --
############################################################################


cdef class Set_generic(Parent): # Cannot use Parent because Element._parent is Parent
    """
    Abstract base class for sets.
    """
    def category(self):
        """
        The category that this set belongs to, which is the category
        of all sets.

        EXAMPLES:
            sage: Set(QQ).category()
            Category of sets
        """
        import sage.categories.all
        return sage.categories.all.Sets()

    def object(self):
        return self

    def __nonzero__(self):
        """
        A set is considered True unless it is empty, in which case it is
        considered to be False.

        EXAMPLES:
            sage: bool(Set(QQ))
            True
            sage: bool(Set(GF(3)))
            True
        """
        return not (self.is_finite() and len(self) == 0)


import types
cdef _type_set_cache = {}

def Set_PythonType(theType):
    """
    Return the (unique) Parent that represents the set of Python objects
    of a specified type.

    EXAMPLES:
        sage: from sage.structure.parent import Set_PythonType
        sage: Set_PythonType(list)
        Set of Python objects of type 'list'
        sage: Set_PythonType(list) is Set_PythonType(list)
        True
        sage: S = Set_PythonType(tuple)
        sage: S([1,2,3])
        (1, 2, 3)
    """
    try:
        return _type_set_cache[theType]
    except KeyError:
        _type_set_cache[theType] = theSet = Set_PythonType_class(theType)
        return theSet

cdef class Set_PythonType_class(Set_generic):

    cdef _type

    def __init__(self, theType):
        from sage.categories.category_types import Sets
        Set_generic.__init__(self, element_constructor=theType, categories=[Sets()])
        self._type = theType

    def __call__(self, x):
        """
        This doesn't return Elements, but actual objects of the given type.

        EXAMPLES:
            sage: S = sage.structure.parent.Set_PythonType(float)
            sage: S(5)
            5.0
            sage: S(9/3)
            3.0
            sage: S(1/3)
            0.333333333333333...
        """
        return self._type(x)

    def __hash__(self):
        """
        TESTS:
            sage: S = sage.structure.parent.Set_PythonType(int)
            sage: hash(S) == -hash(int)
            True
        """
        return -hash(self._type)

    cpdef int _cmp_(self, other) except -100:
        """
        Two Python type sets are considered the same if they contain the same
        type.

        EXAMPLES:
            sage: S = sage.structure.parent.Set_PythonType(int)
            sage: S == S
            True
            sage: S == sage.structure.parent.Set_PythonType(float)
            False
        """
        if self is other:
            return 0
        if isinstance(other, Set_PythonType_class):
            return cmp(self._type, other._type)
        else:
            return cmp(self._type, other)

    def __contains__(self, x):
        """
        Only things of the right type (or subtypes thereof) are considered to
        belong to the set.

        EXAMPLES:
            sage: S = sage.structure.parent.Set_PythonType(tuple)
            sage: (1,2,3) in S
            True
            sage: () in S
            True
            sage: [1,2] in S
            False
        """
        return isinstance(x, self._type)

    def _repr_(self):
        """
        EXAMPLES:
            sage: sage.structure.parent.Set_PythonType(tuple)
            Set of Python objects of type 'tuple'
            sage: sage.structure.parent.Set_PythonType(Integer)
            Set of Python objects of type 'sage.rings.integer.Integer'
            sage: sage.structure.parent.Set_PythonType(Parent)
            Set of Python objects of type 'sage.structure.parent.Parent'
        """
        return "Set of Python objects of %s"%(str(self._type)[1:-1])

    def object(self):
        """
        EXAMPLES:
            sage: S = sage.structure.parent.Set_PythonType(tuple)
            sage: S.object()
            <type 'tuple'>
        """
        return self._type

    def cardinality(self):
        """
        EXAMPLES:
            sage: S = sage.structure.parent.Set_PythonType(bool)
            sage: S.cardinality()
            2
            sage: S = sage.structure.parent.Set_PythonType(int)
            sage: S.cardinality()
            4294967296                        # 32-bit
            18446744073709551616              # 64-bit
            sage: S = sage.structure.parent.Set_PythonType(float)
            sage: S.cardinality()
            18437736874454810627
            sage: S = sage.structure.parent.Set_PythonType(long)
            sage: S.cardinality()
            +Infinity
        """
        from sage.rings.integer import Integer
        two = Integer(2)
        if self._type is bool:
            return two
        elif self._type is int:
            import sys
            return two*sys.maxint + 2
        elif self._type is float:
            return 2 * two**52 * (two**11 - 1) + 3 # all NaN's are the same from Python's point of view
        else:
            # probably
            import sage.rings.infinity
            return sage.rings.infinity.infinity

# These functions are to guarantee that user defined _lmul_, _rmul_,
# _l_action_, _r_action_ do not in turn call __mul__ on their
# arguments, leading to an infinite loop.

cdef object _coerce_test_list = []

class EltPair:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __eq__(self, other):
        return type(self.x) is type(other.x) and self.x is other.x and type(self.y) is type(other.y) and self.y is other.y
    def __repr__(self):
        return "%r (%r), %r (%r)" % (self.x, type(self.x), self.y, type(self.y))

cdef bint _register_pair(x, y) except -1:
    both = EltPair(x,y)
#    print _coerce_test_list, " + ", both
    if both in _coerce_test_list:
#        print "Uh oh..."
#        print _coerce_test_list
#        print both
        xp = type(x) if PY_TYPE_CHECK(x, Parent) else parent_c(x)
        yp = type(y) if PY_TYPE_CHECK(y, Parent) else parent_c(y)
        #print xp, yp
#        if str(y) != "Complex Double Field":
#            raise Exception, "Infinite loop in action of %s (parent %s) and %s (parent %s)!" % (x, xp, y, yp)
        raise NotImplementedError, "Infinite loop in action of %s (parent %s) and %s (parent %s)!" % (x, xp, y, yp)
    _coerce_test_list.append(both)
    return 0

cdef bint _unregister_pair(x, y) except -1:
    try:
        _coerce_test_list.remove(EltPair(x,y))
    except (ValueError, NotImplementedError):
        pass


empty_set = Set_generic()

def normalize_names(ngens, names):
    """
    TESTS:
        sage: sage.structure.parent.normalize_names(5, 'x')
        ('x0', 'x1', 'x2', 'x3', 'x4')
        sage: sage.structure.parent.normalize_names(2, ['x','y'])
        ('x', 'y')
    """
    return empty_set.normalize_names(ngens, names)
