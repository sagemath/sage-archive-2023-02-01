cdef extern from 'symmetrica/def.h':
    INT mult_schubert_schubert(OP a, OP b, OP result)
    INT m_perm_sch(OP a, OP b)
    INT t_SCHUBERT_POLYNOM(OP a, OP b)
    INT t_POLYNOM_SCHUBERT(OP a, OP b)
    INT mult_schubert_variable(OP a, OP i, OP r)
    INT divdiff_perm_schubert(OP perm, OP sb, OP res)
    INT scalarproduct_schubert(OP a, OP b, OP c)
    INT divdiff_schubert(OP a, OP schub, OP res)

    INT t_2SCHUBERT_POLYNOM(OP a,OP b)
    INT mult_schubert_polynom(OP a,OP b,OP c)






def mult_schubert_schubert_symmetrica(a, b):
    """
    Multiplies the Schubert polynomials a and b.
    """
    late_import()

    cdef OP ca = callocobject(), cb = callocobject(), cres = callocobject()

    if isinstance(a, (Permutation_class, builtinlist)) and isinstance(b, (Permutation_class, builtinlist)):
        _op_schubert_perm(a, ca)
        _op_schubert_perm(b, cb)
    else:
        ab = a.parent().base_ring()
        bb = b.parent().base_ring()
        if ab == bb and (ab == QQ or ab == ZZ):
            _op_schubert_sp(a, ca)
            _op_schubert_sp(b, cb)
        else:
            raise ValueError, "a and b must be Schubert polynomials over ZZ or QQ"

    _sig_on
    mult_schubert_schubert(ca, cb, cres)
    _sig_off

    res = _py(cres)

    freeall(ca)
    freeall(cb)
    freeall(cres)

    return res

def t_SCHUBERT_POLYNOM_symmetrica(a):
    late_import()

    cdef OP ca = callocobject(), cres = callocobject()

    if isinstance(a, (Permutation_class, builtinlist)):
        _op_schubert_perm(a, ca)
    else:
        ab = a.parent().base_ring()
        if (ab == QQ or ab == ZZ):
            _op_schubert_sp(a, ca)
        else:
            raise ValueError, "a and b must be Schubert polynomials over ZZ or QQ"


    _sig_on
    t_SCHUBERT_POLYNOM(ca, cres)
    _sig_off

    res = _py(cres)

    freeall(ca)
    freeall(cres)

    return res

def t_POLYNOM_SCHUBERT_symmetrica(a):
    cdef OP ca = callocobject(), cres = callocobject()

    _op_polynom(a, ca)

    _sig_on
    t_POLYNOM_SCHUBERT(ca, cres)
    _sig_off

    res = _py(cres)

    freeall(ca)
    freeall(cres)

    return res

def mult_schubert_variable_symmetrica(a, i):
    late_import()

    cdef OP ca = callocobject(), ci = callocobject(),  cres = callocobject()

    if isinstance(a, (Permutation_class, builtinlist)):
        _op_schubert_perm(a, ca)
    else:
        ab = a.parent().base_ring()
        if (ab == QQ or ab == ZZ):
            _op_schubert_sp(a, ca)
        else:
            raise ValueError, "a and b must be Schubert polynomials over ZZ or QQ"
    _op_integer(i, ci)

    _sig_on
    mult_schubert_variable(ca, ci, cres)
    _sig_off

    res = _py(cres)

    freeall(ca)
    freeall(ci)
    freeall(cres)

    return res


def divdiff_perm_schubert_symmetrica(perm, a):
    late_import()

    cdef OP ca = callocobject(), cperm = callocobject(),  cres = callocobject()

    if isinstance(a, (Permutation_class, builtinlist)):
        _op_schubert_perm(a, ca)
    else:
        ab = a.parent().base_ring()
        if (ab == QQ or ab == ZZ):
            _op_schubert_sp(a, ca)
        else:
            raise ValueError, "a and b must be Schubert polynomials over ZZ or QQ"
    _op_permutation(perm, cperm)

    _sig_on
    divdiff_perm_schubert(cperm, ca, cres)
    _sig_off

    res = _py(cres)

    freeall(ca)
    freeall(cperm)
    freeall(cres)

    return res


def scalarproduct_schubert_symmetrica(a, b):
    late_import()

    cdef OP ca = callocobject(), cb = callocobject(), cres = callocobject()

    if isinstance(a, (Permutation_class, builtinlist)) and isinstance(b, (Permutation_class, builtinlist)):
        _op_schubert_perm(a, ca)
        _op_schubert_perm(b, cb)
    else:
        ab = a.parent().base_ring()
        bb = b.parent().base_ring()
        if ab == bb and (ab == QQ or ab == ZZ):
            _op_schubert_sp(a, ca)
            _op_schubert_sp(b, cb)
        else:
            raise ValueError, "a and b must be Schubert polynomials over ZZ or QQ"


    _sig_on
    scalarproduct_schubert(ca, cb, cres)
    _sig_off


    if empty_listp(cres):
        freeall(ca)
        freeall(cb)
        freeall(cres)
        return Integer(0)

    res = _py(cres)

    freeall(ca)
    freeall(cb)
    freeall(cres)

    return res

def divdiff_schubert_symmetrica(i, a):
    late_import()

    cdef OP ca = callocobject(), ci = callocobject(),  cres = callocobject()

    if isinstance(a, (Permutation_class, builtinlist)):
        _op_schubert_perm(a, ca)
    else:
        ab = a.parent().base_ring()
        if (ab == QQ or ab == ZZ):
            _op_schubert_sp(a, ca)
        else:
            raise ValueError, "a and b must be Schubert polynomials over ZZ or QQ"
    _op_integer(i, ci)

    _sig_on
    divdiff_schubert(ci, ca, cres)
    _sig_off

    res = _py(cres)

    freeall(ca)
    freeall(ci)
    freeall(cres)

    return res
