from sage.libs.ntl.types cimport ZZ_c

cdef extern from "eclib/interface.h":
    ctypedef struct bigint:  #eclib uses NTL in Sage -- we call Cremona's "bigint" ZZ_c.
        pass
    ZZ_c new_bigint "to_ZZ"(long)
    int I2int(ZZ_c)

cdef extern from "eclib/bigrat.h":
    ctypedef struct bigrational "bigrational":
        pass
    bigrational new_bigrational "bigrational"(ZZ_c num, ZZ_c den)
    cdef ZZ_c bigrational_num "num"(bigrational q)
    cdef ZZ_c bigrational_den "den"(bigrational q)

cdef extern from "eclib/rat.h":
    ctypedef struct rational "rational":
        pass
    rational new_rational "rational"(long num, long den)
    cdef long rational_num "num"(rational q)
    cdef long rational_den "den"(rational q)

cdef extern from "eclib/curve.h":
    ctypedef struct Curve:
        void (* getai)(ZZ_c a1, ZZ_c a2, ZZ_c a3, ZZ_c a4, ZZ_c a6)
    Curve *new_Curve "new Curve"( ZZ_c aa1,  ZZ_c aa2,  ZZ_c aa3,  ZZ_c aa4,  ZZ_c aa6)

    ctypedef struct Curvedata:
        void (* getai)(ZZ_c a1, ZZ_c a2, ZZ_c a3, ZZ_c a4, ZZ_c a6)
    Curvedata *new_Curvedata "new Curvedata"( Curve C, int m)

    ctypedef struct CurveRed:
        pass
    CurveRed *new_CurveRed "new CurveRed"( Curvedata CD)

    ZZ_c getconductor(CurveRed CR)


