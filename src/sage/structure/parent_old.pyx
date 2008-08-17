r"""
Base class for parent objects

CLASS HIEARCHY:

SageObject
    Parent
        ParentWithBase
            ParentWithGens


TESTS:
This came up in some subtle bug once.
    sage: gp(2) + gap(3)
    5
"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport sage_object
import operator
from parent import Set_PythonType, Set_PythonType_class

include '../ext/python_object.pxi'
include '../ext/python_bool.pxi'
include '../ext/stdsage.pxi'

def is_Parent(x):
    """
    Return True if x is a parent object, i.e., derives from
    sage.structure.parent.Parent and False otherwise.

    EXAMPLES:
        sage: is_Parent(2/3)
        False
        sage: is_Parent(ZZ)
        True
        sage: is_Parent(Primes())
        True
    """
    return PyBool_FromLong(PyObject_TypeCheck(x, Parent))

cdef inline check_old_coerce(Parent p):
    if p._element_constructor is not None:
        raise RuntimeError, "%s still using old coercion framework" % p


## def make_parent_v0(_class, _dict, has_coerce_map_from):
##     """
##     This should work for any Python class deriving from this, as long
##     as it doesn't implement some screwy __new__() method.
##     """
##     cdef Parent new_object
##     new_object = _class.__new__(_class)
##     if not _dict is None:
##         new_object.__dict__ = _dict
##     new_object._has_coerce_map_from = has_coerce_map_from
##     return new_object

cdef class Parent(parent.Parent):
    """
    Parents are the SAGE/mathematical analogues of container objects
    in computer science.
    """

    def __init__(self, coerce_from=[], actions=[], embeddings=[]):
        # TODO: many classes don't call this at all, but __new__ crashes SAGE
#        if len(coerce_from) > 0:
#            print type(self), coerce_from
        self._coerce_from_list = list(coerce_from)
        self._coerce_from_hash = {}
        self._action_list = list(actions)
        self._action_hash = {}

        cdef parent.Parent other
        for mor in embeddings:
            other = mor.domain()
            print "embedding", self, " --> ", other
            print mor
            other.init_coerce() # TODO remove when we can
            other._coerce_from_list.append(mor)

        # old
        self._has_coerce_map_from = {}

    cdef int init_coerce(self, bint warn=False) except -1:
        parent.Parent.init_coerce(self, warn)


    #################################################################################
    # New Coercion support functionality
    #################################################################################

#    def coerce_map_from(self, S):
#        return self.coerce_map_from_c(S)

    cpdef coerce_map_from_c(self, S):
        check_old_coerce(self)
        if S is self:
            from sage.categories.homset import Hom
            return Hom(self, self).identity()
        elif S == self:
            # non-unique parents
            from sage.categories.homset import Hom
            from sage.categories.morphism import CallMorphism
            return CallMorphism(Hom(S, self))
        elif isinstance(S, Set_PythonType_class):
            return self.coerce_map_from_c(S._type)
        if self._coerce_from_hash is None: # this is because parent.__init__() does not always get called
            self.init_coerce()
        cdef object ret
        try:
            ret = PyObject_GetItem(self._coerce_from_hash,S)
            return ret
        except KeyError:
            pass
        if HAS_DICTIONARY(self):
            mor = self.coerce_map_from_impl(S)
        else:
            mor = self.coerce_map_from_c_impl(S)
        import sage.categories.morphism
        import sage.categories.map
        if mor is True:
            mor = sage.categories.morphism.CallMorphism(S, self)
        elif mor is False:
            mor = None
        elif mor is not None and not isinstance(mor, sage.categories.map.Map):
            raise TypeError, "coerce_map_from_impl must return a boolean, None, or an explicit Map"
        if mor is not None:
            self._coerce_from_hash[S] = mor # TODO: if this is None, could it be non-None in the future?
        return mor

    def coerce_map_from_impl(self, S):
        check_old_coerce(self)
        return self.coerce_map_from_c_impl(S)

    cdef coerce_map_from_c_impl(self, S):
        check_old_coerce(self)
        import sage.categories.morphism
        from sage.categories.map import Map
        from sage.categories.homset import Hom
        cdef parent.Parent R
        for mor in self._coerce_from_list:
            if PY_TYPE_CHECK(mor, Map):
                R = mor.domain()
            else:
                R = mor
                mor = sage.categories.morphism.CallMorphism(Hom(R, self))
                i = self._coerce_from_list.index(R)
                self._coerce_from_list[i] = mor # cache in case we need it again
            if R is S:
                return mor
            else:
                connecting = R.coerce_map_from(S)
                if connecting is not None:
                    return mor * connecting

        # Piggyback off the old code for now
        # WARNING: when working on this, make sure circular dependancies aren't introduced!
        if self.has_coerce_map_from_c(S):
            if isinstance(S, type):
                S = Set_PythonType(S)
            return sage.categories.morphism.CallMorphism(Hom(S, self))
        else:
            return None

#    def get_action(self, S, op=operator.mul, self_on_left=True):
#        return self.get_action_c(S, op, self_on_left)

    cpdef get_action_c(self, S, op, bint self_on_left):
        check_old_coerce(self)
        try:
            if self._action_hash is None: # this is because parent.__init__() does not always get called
                self.init_coerce()
            return self._action_hash[S, op, self_on_left]
        except KeyError:
            pass
        if HAS_DICTIONARY(self):
            action = self.get_action_impl(S, op, self_on_left)
        else:
            action = self.get_action_c_impl(S, op, self_on_left)
        if action is not None:
            from sage.categories.action import Action
            if not isinstance(action, Action):
                raise TypeError, "get_action_impl must return None or an Action"
            self._action_hash[S, op, self_on_left] = action
        return action

    def get_action_impl(self, S, op, self_on_left):
        check_old_coerce(self)
        return self.get_action_c_impl(S, op, self_on_left)

    cdef get_action_c_impl(self, S, op, bint self_on_left):
        check_old_coerce(self)
        # G acts on S, G -> G', R -> S => G' acts on R (?)
        from sage.categories.action import Action, PrecomposedAction
        from sage.categories.homset import Hom
        from coerce_actions import LeftModuleAction, RightModuleAction
        cdef parent.Parent R
        for action in self._action_list:
            if PY_TYPE_CHECK(action, Action):
                if self_on_left:
                    if action.left() is not self: continue
                    R = action.right()
                else:
                    if action.right() is not self: continue
                    R = action.left()
            elif op is operator.mul:
                try:
                    R = action
                    _register_pair(x,y) # to kill circular recursion
                    if self_on_left:
                        action = LeftModuleAction(S, self) # self is acted on from right
                    else:
                        action = RightModuleAction(S, self) # self is acted on from left
                    _unregister_pair(x,y)
                    i = self._action_list.index(R)
                    self._action_list[i] = action
                except TypeError:
                    continue
            else:
                continue # only try mul if not specified
            if R is S:
                return action
            else:
                connecting = R.coerce_map_from_c(S) # S -> R
                if connecting is not None:
                    if self_on_left:
                        return PrecomposedAction(action, None, connecting)
                    else:
                        return PrecomposedAction(action, connecting, None)


        if op is operator.mul and PY_TYPE_CHECK(S, parent.Parent):
            from coerce_actions import LeftModuleAction, RightModuleAction, LAction, RAction
            # Actors define _l_action_ and _r_action_
            # Acted-on elements define _lmul_ and _rmul_

            # TODO: if _xmul_/_x_action_ code does stuff like
            # if self == 0:
            #    return self
            # then _an_element_c() == 0 could be very bad.
            #
            #
            x = self.an_element()
            y = (<parent.Parent>S).an_element()
#            print "looking for action ", self, "<--->", S
#            print "looking for action ", x, "<--->", y

            _register_pair(x,y) # this is to avoid possible infinite loops
            if self_on_left:
                try:
#                    print "RightModuleAction"
                    action = RightModuleAction(S, self) # this will test _lmul_
                    _unregister_pair(x,y)
#                    print "got", action
                    return action
                except (NotImplementedError, TypeError, AttributeError, ValueError), ex:
                    pass

                try:
#                    print "LAction"
                    z = x._l_action_(y)
                    _unregister_pair(x,y)
                    return LAction(self, S)
                except (NotImplementedError, TypeError, AttributeError, ValueError):
                    pass

            else:
                try:
#                    print "LeftModuleAction"
                    action = LeftModuleAction(S, self) # this will test _rmul_
                    _unregister_pair(x,y)
#                    print "got", action
                    return action
                except (NotImplementedError, TypeError, AttributeError, ValueError), ex:
                    pass

                try:
#                    print "RAction"
                    z = x._r_action_(y)
                    _unregister_pair(x,y)
                    return RAction(self, S)
                except (NotImplementedError, TypeError, AttributeError, ValueError):
                    pass


            try:
                # maybe there is a more clever way of detecting ZZ than importing here...
                from sage.rings.integer_ring import ZZ
                if S is ZZ and not self.has_coerce_map_from(ZZ):
#                    print "IntegerMulAction"
                    from sage.structure.coerce import IntegerMulAction
                    action = IntegerMulAction(S, self, not self_on_left)
#                    print "got", action
                    _unregister_pair(x,y)
                    return action
            except (NotImplementedError, TypeError, AttributeError, ValueError):
                pass

            _unregister_pair(x,y)

#            print "found nothing"


    #################################################################################
    # Coercion support functionality
    #################################################################################

    def _coerce_(self, x):            # Call this from Python (do not override!)
        check_old_coerce(self)
        return self._coerce_c(x)

    cpdef _coerce_c(self, x):          # DO NOT OVERRIDE THIS (call it)
        check_old_coerce(self)
        try:
            P = x.parent()   # todo -- optimize
            if P is self:
                return x
        except AttributeError, msg:
            pass
        if HAS_DICTIONARY(self):
            return self._coerce_impl(x)
        else:
            return self._coerce_c_impl(x)

    cdef _coerce_c_impl(self, x):     # OVERRIDE THIS FOR SAGEX CLASES
        """
        Canonically coerce x in assuming that the parent of x is not
        equal to self.
        """
        check_old_coerce(self)
        raise TypeError

    def _coerce_impl(self, x):        # OVERRIDE THIS FOR PYTHON CLASSES
        """
        Canonically coerce x in assuming that the parent of x is not
        equal to self.
        """
        check_old_coerce(self)
        return self._coerce_c_impl(x)

    def _coerce_try(self, x, v):
        """
        Given a list v of rings, try to coerce x canonically into each
        one in turn.  Return the __call__ coercion of the result into
        self of the first canonical coercion that succeeds.  Raise a
        TypeError if none of them succeed.

        INPUT:
             x -- Python object
             v -- parent object or list (iterator) of parent objects
        """
        check_old_coerce(self)
        if not isinstance(v, list):
            v = [v]

        for R in v:
            try:
                y = R._coerce_(x)
                return self(y)
            except (TypeError, AttributeError), msg:
                pass
        raise TypeError, "no canonical coercion of element into self"

    def _coerce_self(self, x):
        check_old_coerce(self)
        return self._coerce_self_c(x)

    cdef _coerce_self_c(self, x):
        """
        Try to canonically coerce x into self.
        Return result on success or raise TypeError on failure.
        """
        check_old_coerce(self)
        # todo -- optimize?
        try:
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return self(x)
        except AttributeError:
            pass
        raise TypeError, "no canonical coercion to self defined"

#    def has_coerce_map_from(self, S):
#        return self.has_coerce_map_from_c(S)

    cpdef has_coerce_map_from_c(self, S):
        """
        Return True if there is a natural map from S to self.
        Otherwise, return False.
        """
        check_old_coerce(self)
        if self == S:
            return True
        if self._has_coerce_map_from is None:
            self._has_coerce_map_from = {}
        else:
            try:
                return self._has_coerce_map_from[S]
            except KeyError:
                pass
        if HAS_DICTIONARY(self):
            x = self.has_coerce_map_from_impl(S)
        else:
            x = self.has_coerce_map_from_c_impl(S)
        self._has_coerce_map_from[S] = x
        return x

    def has_coerce_map_from_impl(self, S):
        check_old_coerce(self)
        return self.has_coerce_map_from_c_impl(S)

    cdef has_coerce_map_from_c_impl(self, S):
        check_old_coerce(self)
        if not PY_TYPE_CHECK(S, parent.Parent):
            return False
        try:
            self._coerce_c((<parent.Parent>S).an_element())
        except TypeError:
            return False
        except NotImplementedError, msg:
            raise NotImplementedError, "%s\nAlso, please make sure you have implemented has_coerce_map_from_impl or has_coerce_map_from_c_impl (or better _an_element_c_impl or _an_element_impl if possible) for %s"%(msg,self)
        return True

    def _an_element_impl(self):     # override this in Python
        check_old_coerce(self)
        return self._an_element_c_impl()

    cpdef _an_element_c_impl(self):  # override this in SageX
        """
        Returns an element of self. Want it in sufficent generality
        that poorly-written functions won't work when they're not
        supposed to. This is cached so doesn't have to be super fast.
        """
        check_old_coerce(self)
        try:
            return self.gen(0)
        except:
            pass

        try:
            return self.gen()
        except:
            pass

        from sage.rings.infinity import infinity
        for x in ['_an_element_', 'pi', 1.2, 2, 1, 0, infinity]:
            try:
                return self(x)
            except (TypeError, NameError, NotImplementedError, AttributeError, ValueError):
                pass

        raise NotImplementedError, "please implement _an_element_c_impl or _an_element_impl for %s"%self

    def _an_element(self):        # do not override this (call from Python)
        check_old_coerce(self)
        return self._an_element_c()

    cpdef _an_element_c(self):     # do not override this (call from SageX)
        check_old_coerce(self)
        if not self.__an_element is None:
            return self.__an_element
        if HAS_DICTIONARY(self):
            self.__an_element = self._an_element_impl()
        else:
            self.__an_element = self._an_element_c_impl()
        return self.__an_element


    ################################################
    # Comparison of parent objects
    ################################################
    cdef _richcmp(left, right, int op):
        """
        Compare left and right.
        """
        check_old_coerce(left)
        cdef int r

        if not PY_TYPE_CHECK(right, parent.Parent) or not PY_TYPE_CHECK(left, parent.Parent):
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
                r = left._cmp_c_impl(right)

        if op == 0:  #<
            return PyBool_FromLong(r  < 0)
        elif op == 2: #==
            return PyBool_FromLong(r == 0)
        elif op == 4: #>
            return PyBool_FromLong(r  > 0)
        elif op == 1: #<=
            return PyBool_FromLong(r <= 0)
        elif op == 3: #!=
            return PyBool_FromLong(r != 0)
        elif op == 5: #>=
            return PyBool_FromLong(r >= 0)

##     ####################################################################
##     # For a derived SageX class, you **must** put the following in
##     # your subclasses, in order for it to take advantage of the
##     # above generic comparison code.  You must also define
##     # _cmp_c_impl for a SageX class.
##     #
##     # For a derived Python class, simply define __cmp__.
##     ####################################################################
##     def __richcmp__(left, right, int op):
##         return (<Parent>left)._richcmp(right, op)

##         # NOT NEEDED, since all attributes are public!
##     def __reduce__(self):
##         if HAS_DICTIONARY(self):
##             _dict = self.__dict__
##         else:
##             _dict = None
##         return (make_parent_v0, (self.__class__, _dict, self._has_coerce_map_from))

    cdef int _cmp_c_impl(left, parent.Parent right) except -2:
        check_old_coerce(left)
        pass
        # this would be nice to do, but we can't since
        # it leads to infinite recurssions -- and is slow -- and this
        # stuff must be fast!
        #if right.has_coerce_map_from(left):
        #    if left.has_coerce_map_from(right):
        #        return 0
        #    else:
        #        return -1
        if (<PyObject*>left) < (<PyObject*>right):
            return -1
        elif (<PyObject*>left) > (<PyObject*>right):
            return 1
        return 0

##     def __cmp__(left, right):
##         return left._cmp_c_impl(right)   # default


    ############################################################################
    # Coercion Compatibility Layer
    ############################################################################

    cpdef _coerce_map_from_(self, S):
        if self._element_constructor is None:
            return self.coerce_map_from_c(S)
        else:
            return None

    cpdef _get_action_(self, other, op, bint self_on_left):
        if self._element_constructor is None:
            return self.get_action_c(other, op, self_on_left)
        else:
            return None

    cpdef _an_element_(self):
        if self._element_constructor is None:
            return self._an_element_c()
        else:
            return parent.Parent._an_element_(self)

    cpdef _generic_convert_map(self, S):
        if self._element_constructor is None:
            from sage.categories.morphism import CallMorphism
            from sage.categories.homset import Hom
            return CallMorphism(Hom(S, self))
        else:
            return parent.Parent._generic_convert_map(self, S)


############################################################################
# Set baseclass --
############################################################################

# These functions are to guerentee that user defined _lmul_, _rmul_, _l_action_, _r_action_ do
# not in turn call __mul__ on their arguments, leading to an infinite loop.

cdef object _coerce_test_list = []

class EltPair:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __eq__(self, other):
        return type(self.x) is type(other.x) and self.x == other.x and type(self.y) is type(other.y) and self.y == other.y
    def __repr__(self):
        return "%r (%r), %r (%r)" % (self.x, type(self.x), self.y, type(self.y))

cdef bint _register_pair(x, y) except -1:
    both = EltPair(x,y)
#    print _coerce_test_list, " + ", both
    if both in _coerce_test_list:
#        print "Uh oh..."
#        print _coerce_test_list
#        print both
        raise NotImplementedError, "Infinite loop in multiplication of %s (parent %s) and %s (parent %s)!" % (x, x.parent(), y, y.parent())
    _coerce_test_list.append(both)
    return 0

cdef void _unregister_pair(x, y):
    try:
        _coerce_test_list.remove(EltPair(x,y))
    except (ValueError, NotImplementedError):
        pass
