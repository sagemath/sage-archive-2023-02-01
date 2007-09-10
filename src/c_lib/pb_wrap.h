#include "polybori/polybori.h"
#include "ccobject.h"
#include <sstream>

USING_NAMESPACE_PBORI


BoolePolyRing* BPRing_Construct(void* mem, BoolePolyRing::size_type nvars,
        BoolePolyRing::ordercode_type order){
    return new(mem) BoolePolyRing(nvars, order);
}

BoolePolynomial* BPolyNewDD(const BoolePolyRing::dd_type &d){
    return new BoolePolynomial(d);
}

BoolePolynomial* BPolyNewBPoly(const BoolePolynomial &d){
    return new BoolePolynomial(d);
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

BoolePolynomial* PBPoly_add(BoolePolynomial* rval, BoolePolynomial* lval)
{
    BoolePolynomial* n = new BoolePolynomial(*rval + *lval);
    return n;
}

BoolePolynomial* PBPoly_mul(BoolePolynomial* rval, BoolePolynomial* lval)
{
    BoolePolynomial* n = new BoolePolynomial(*rval * *lval);
    return n;
}
