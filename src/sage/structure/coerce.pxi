# These are shared inline functions.

#################################################################################
# fast tests to verify no coersion is needed
#################################################################################

cdef inline bint have_same_parent(left, right):
    """
    Return nonzero true value if and only if left and right are
    elements and have the same parent.
    """
    # (We know at least one of the arguments is an Element. So if
    # their types are *equal* (fast to check) then they are both
    # Elements.  Otherwise use the slower test via PY_TYPE_CHECK.)
    if PY_TYPE(left) is PY_TYPE(right):
        return (<Element>left)._parent is (<Element>right)._parent

    if PY_TYPE_CHECK(right, Element) and PY_TYPE_CHECK(left, Element):
        return (<Element>left)._parent is (<Element>right)._parent

    return 0

cdef inline _verify_canonical_coercion_c(x, y):
    if not have_same_parent(x,y):
        raise RuntimeError, """There is a bug in the coercion code in SAGE.
Both x (=%s) and y (=%s) are supposed to have identical parents but they don't.
In fact, x has parent '%s'
whereas y has parent '%s'"""%(x,y,parent_c(x),parent_c(y))
    return x, y

cdef inline bint have_same_base(Element x, Element y):
    return x._parent._base is y._parent._base

#################################################################################
# parent
#################################################################################
cdef inline parent_c(x):
    if PY_TYPE_CHECK(x,Element):
        return (<Element>x)._parent
    elif hasattr(x, 'parent'):
        return x.parent()
    else:
        return <object>PY_TYPE(x)

def parent(x):
    return parent_c(x)

#################################################################################
# errors
#################################################################################
cdef D = {'mul':'*', 'add':'+', 'sub':'-', 'div':'/'}

cdef inline arith_error_message(x, y, op):
        try:
            n = D[op.__name__]
        except KeyError:
            n = op.__name__
        return "unsupported operand parent(s) for '%s': '%s' and '%s'"%(n, parent_c(x), parent_c(y))

cdef inline ModuleElement _add_c(ModuleElement left, ModuleElement right):
    if HAS_DICTIONARY(left):
        return left._add_(right)
    else:
        return left._add_c_impl(right)

cdef inline ModuleElement _iadd_c(ModuleElement left, ModuleElement right):
    if HAS_DICTIONARY(left):
        return left._iadd_(right)
    else:
        return left._iadd_c_impl(right)
