import sage.libs.pari.gen as gen
import sage.misc.misc

from math import log
from sage.rings.all import *
I = ComplexField().gen()

def pari(x):
    return gen.pari(x)

def python(z):
    t = z.type()
    if t == "t_REAL":
        wordsize = 32
        if sage.misc.misc.is_64_bit:
            wordsize = 64
        prec = (z.gen_length() - 2) * wordsize
        R = RealField(prec)
        return R(z)
    elif t == "t_FRAC":
         Q = RationalField()
         return Q(z)
    elif t == "t_COMPLEX":
        # The components of a Pari COMPLEX might not be floats.
        # We compute the smaller of the real and imaginary precisions,
        # treating non-floats as infinitely precise and defaulting
        # to the default 53-bit precision if both components are
        # non-floats.
        wordsize = 32
        if sage.misc.misc.is_64_bit:
            wordsize = 64
        prec = None
        re = z.real()
        if re.type() == "t_REAL":
            prec = (re.gen_length() - 2) * wordsize
        im = z.imag()
        if im.type() == "t_REAL":
            im_prec = (im.gen_length() - 2) * wordsize
            if prec is None or im_prec < prec:
                prec = im_prec
        if prec is None:
            fld = ComplexField()
        else:
            fld = ComplexField(prec)
        return fld(re, im)
    elif t == "t_VEC":
        return [python(x) for x in z.python_list()]
    elif t == "t_VECSMALL":
        return [ZZ(x) for x in z.python_list_small()]
    elif t == "t_MAT":
        return [python(z[i,j]) for i in range(z.nrows()) for j in range(z.ncols())]
    else:
        return eval(str(z))
