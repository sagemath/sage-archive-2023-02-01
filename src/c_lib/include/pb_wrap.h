#include "polybori.h"
#include "groebner_alg.h"
#include "nf.h"
#include "ccobject.h"
#include <sstream>

USING_NAMESPACE_PBORI
USING_NAMESPACE_PBORIGB


BoolePolyRing* PBRing_construct(void* mem, BoolePolyRing::size_type nvars,
        BoolePolyRing::ordercode_type order){
    return new(mem) BoolePolyRing(nvars, order);
}

CCuddNavigator* PBNavigator_construct(void* mem, const CCuddNavigator &d){
    return new(mem) CCuddNavigator(d);
}

BooleSet* PBSet_construct_dd(void* mem, const BooleSet::dd_type &d) {
    return new(mem) BooleSet(d);
}

BooleSet* PBSet_construct_pbnav(void* mem, const CCuddNavigator &d) {
    return new(mem) BooleSet(d);
}

BoolePolynomial* PBPoly_construct_dd(void* mem, const BoolePolyRing::dd_type &d) {
    return new(mem) BoolePolynomial(d);
}

BoolePolynomial* PBPoly_construct_pbset(void* mem, const BooleSet &d) {
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

BooleMonomial* PBMonom_construct_pbmonom(void *mem, const BooleMonomial &d) {
    return new(mem) BooleMonomial(d);
}

BoolePolynomial* PBMonom_construct_dd(void* mem, const BooleMonomial::dd_type &d) {
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
