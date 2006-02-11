import sage.libs.pari.gen as gen

from math import log
from sage.rings.all import *
I = ComplexField().gen()

def pari(x):
    return gen.pari(x)

def python(z):
    t = z.type()
    if t == "t_REAL":
        s = str(z).replace(" E","e")
        prec = int(len(s)*3)  # underestimate
        R = RealField(prec)
        return R(s)
    elif t == "t_COMPLEX":
        s = str(z).replace(" E","e")
        return eval(s)
    elif t == "t_VEC":
        return [python(x) for x in z.python_list()]
    elif t == "t_MAT":
        return [python(z[i,j]) for i in range(z.nrows()) for j in range(z.ncols())]
    else:
        return eval(str(z))
