import sage.libs.pari.gen as gen

from math import log
from sage.rings.all import *
I = ComplexField().gen()

def pari(x):
    """
    Return the pari object constructed from a Sage object.

    The work is done by the __call__ method of the class PariInstance,
    which in turn passes the work to any class which has its own
    method _pari_().

    EXAMPLES:
        sage: a = pari(1); a, a.type()
        (1, 't_INT')
        sage: a = pari(1/2); a, a.type()
        (1/2, 't_FRAC')
        sage: a = pari(1/2); a, a.type()
        (1/2, 't_FRAC')

    Conversion from reals uses the real's own precision, here 53 bits (the deault):
        sage: a = pari(1.2); a, a.type(), a.precision()
        (1.20000000000000, 't_REAL', 4) # 32-bit
        (1.20000000000000, 't_REAL', 3) # 64-bit

    Conversion from strings uses the current pari real precision.  By
    default this is 4 words, 38 digits, 128 bits on 64-bit machines
    and 5 words, 19 digits, 64 bits on 32-bit machines.

        sage: a = pari('1.2'); a, a.type(), a.precision()
        (1.20000000000000, 't_REAL', 5) # 32-bit
        (1.20000000000000, 't_REAL', 4) # 64-bit

    Conversion from matrices is supported , but not from vectors; use
    lists instead:

        sage: a = pari(matrix(2,3,[1,2,3,4,5,6])); a, a.type()
        ([1, 2, 3; 4, 5, 6], 't_MAT')

        sage: v = vector([1.2,3.4,5.6])
        sage: v.pari()
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.modules.free_module_element.FreeModuleElement' object has no attribute 'pari'
        sage: b = pari(list(v)); b,b.type()
        ([1.20000000000000, 3.40000000000000, 5.60000000000000], 't_VEC')

    Some more exotic examples:
        sage: K.<a>=NumberField(x^3-2)
        sage: pari(K)
        [x^3 - 2, [1, 1], -108, 1, [[1, 1.25992104989487, 1.58740105196820; 1, -0.629960524947437 - 1.09112363597172*I, -0.793700525984100 + 1.37472963699860*I], [1, 1.25992104989487, 1.58740105196820; 1, -1.72108416091916, 0.581029111014503; 1, 0.461163111024285, -2.16843016298270], 0, [3, 0, 0; 0, 0, 6; 0, 6, 0], [6, 0, 0; 0, 6, 0; 0, 0, 3], [2, 0, 0; 0, 0, 1; 0, 1, 0], [2, [0, 0, 2; 1, 0, 0; 0, 1, 0]]], [1.25992104989487, -0.629960524947437 - 1.09112363597172*I], [1, x, x^2], [1, 0, 0; 0, 1, 0; 0, 0, 1], [1, 0, 0, 0, 0, 2, 0, 2, 0; 0, 1, 0, 1, 0, 0, 0, 0, 2; 0, 0, 1, 0, 1, 0, 1, 0, 0]]

        sage: E = EllipticCurve('37a1')
        sage: pari(E)
        [0, 0, 1, -1, 0, 0, -2, 1, -1, 48, -216, 37, 110592/37, [0.837565435283323, 0.269594436405445, -1.10715987168877]~, 2.99345864623196, 2.45138938198679*I, -0.471319277956811, -1.43545651866868*I, 7.33813274078958]
    """
    return gen.pari(x)

def python(z):
    """
    Return the closest python/Sage equivalent of the given pari object.

    The component parts of a t_COMPLEX may be t_INT, t_REAL, t_INTMOD,
    t_FRAC, t_PADIC.  Th ecomponents need not have the same type
    (e.g. if z=2+1.2*I then z.real() is t_INT while z.imag() is
    t_REAL().  They are converted as follows:

    t_INT:    ZZ[i]
    t_FRAC:   QQ(i)
    t_REAL:   ComplexField(prec) for equivalent precision
    t_INTMOD, t_PADIC: raise NotImplementedError

    EXAMPLES:
        sage: a = pari('(3+I)').python(); a
        i + 3
        sage: a.parent()
        Maximal Order in Number Field in i with defining polynomial x^2 + 1

        sage: a = pari('2^31-1').python(); a
        2147483647
        sage: a.parent()
        Integer Ring

        sage: a = pari('12/34').python(); a
        6/17
        sage: a.parent()
        Rational Field

        sage: a = pari('1.234').python(); a
        1.234000000000000000000000000           # 32-bit
        1.2340000000000000000000000000000000000 # 64-bit
        sage: a.parent()
        Real Field with 96 bits of precision    # 32-bit
        Real Field with 128 bits of precision   # 64-bit

        sage: a = pari('(3+I)/2').python(); a
        1/2*i + 3/2
        sage: a.parent()
        Number Field in i with defining polynomial x^2 + 1

    Conversion of complex numbers: the next example is converting from
    an element of the Symbolic Ring, which goes via the string
    representation and hence the precision is architecture-dependent:
        sage: a = pari(1.0+2.0*I).python(); a
        1.000000000000000000000000000 + 2.000000000000000000000000000*I  # 32-bit
        1.0000000000000000000000000000000000000 + 2.0000000000000000000000000000000000000*I # 64-bit
        sage: type(a)
        <type 'sage.rings.complex_number.ComplexNumber'>
        sage: a.parent()
        Complex Field with 96 bits of precision # 32-bit
        Complex Field with 128 bits of precision # 64-bit

    For architecture-independent complex numbers, start from a
    suitable ComplexField:
        sage: z = pari(CC(1.0+2.0*I)); z
        1.00000000000000 + 2.00000000000000*I
        sage: a=z.python(); a
        1.00000000000000000 + 2.00000000000000000*I
        sage: a.parent()
        Complex Field with 64 bits of precision

        sage: a = pari('[1,2,3,4]')
        sage: a
        [1, 2, 3, 4]
        sage: a.type()
        't_VEC'
        sage: b = a.python(); b
        [1, 2, 3, 4]
        sage: type(b)
        <type 'list'>

        sage: a = pari('[1,2;3,4]')
        sage: a.type()
        't_MAT'
        sage: b = a.python(); b
        [1 2]
        [3 4]
        sage: b.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

    """
    t = z.type()
    if t == "t_REAL":
        return RealField(gen.prec_words_to_bits(z.precision()))(z)
    elif t == "t_FRAC":
         Q = RationalField()
         return Q(z)
    elif t == "t_INT":
         Z = IntegerRing()
         return Z(z)
    elif t == "t_COMPLEX":
        tx = z.real().type()
        ty = z.imag().type()
        if tx in ["t_INTMOD", "t_PADIC"] or ty in ["t_INTMOD", "t_PADIC"]:
            raise NotImplementedError, "No conversion to python available for t_COMPLEX with t_INTMOD or t_PADIC components"
        if tx == "t_REAL" or ty == "t_REAL":
            xprec = z.real().precision() # will be 0 if exact
            yprec = z.imag().precision() # will be 0 if exact
            if xprec == 0:
                prec = gen.prec_words_to_bits(yprec)
            elif yprec == 0:
                prec = gen.prec_words_to_bits(xprec)
            else:
                prec = min(gen.prec_words_to_bits(xprec),gen.prec_words_to_bits(yprec))
            R = RealField(prec)
            C = ComplexField(prec)
            return C(R(z.real()), R(z.imag()))
        if tx == "t_FRAC" or ty == "t_FRAC":
            return QuadraticField(-1,'i')([python(c) for c in list(z)])
        if tx == "t_INT" or ty == "t_INT":
            return QuadraticField(-1,'i').ring_of_integers()([python(c) for c in list(z)])
        raise NotImplementedError, "No conversion to python available for t_COMPLEX with components %s"%(tx,ty)
    elif t == "t_VEC":
        return [python(x) for x in z.python_list()]
    elif t == "t_VECSMALL":
        return [IntegerRing(x) for x in z.python_list_small()]
    elif t == "t_MAT":
        from sage.matrix.constructor import matrix
        return matrix(z.nrows(),z.ncols(),[python(z[i,j]) for i in range(z.nrows()) for j in range(z.ncols())])
    else:
        return eval(str(z))
