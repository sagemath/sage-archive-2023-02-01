#include "polybori.h"
#include "groebner_alg.h"
#include "ccobject.h"
#include <sstream>

USING_NAMESPACE_PBORI
USING_NAMESPACE_PBORIGB


BoolePolyRing* PBRing_construct(void* mem, BoolePolyRing::size_type nvars,
        BoolePolyRing::ordercode_type order){
    return new(mem) BoolePolyRing(nvars, order);
}

BoolePolynomial* PBPoly_construct_dd(void* mem, const BoolePolyRing::dd_type &d) {
    return new(mem) BoolePolynomial(d);
}

BoolePolynomial* PBPoly_construct_pbpoly(void *mem, const BoolePolynomial &d) {
    return new(mem) BoolePolynomial(d);
}

BoolePolynomial* PBPoly_construct_pbmonom(void *mem, const BooleMonomial &d) {
    return new(mem) BoolePolynomial(d);
}

BoolePolynomial* PBPoly_construct_int(void *mem, const int d) {
    return new(mem) BoolePolynomial(d);
}


template<class T>
char* to_str(const T* x)
{
    std::ostringstream instore;
    instore << (*x);
    int n = strlen(instore.str().data());
    char* buf = (char*)malloc(n+1);
    strcpy(buf, instore.str().data());
    return buf;
}

inline BoolePolynomial PBPoly_add(BoolePolynomial rval, BoolePolynomial lval)
{
    return BoolePolynomial(rval + lval);
}

inline BoolePolynomial PBPoly_mul(BoolePolynomial rval, BoolePolynomial lval)
{
    return BoolePolynomial(rval * lval);
}
