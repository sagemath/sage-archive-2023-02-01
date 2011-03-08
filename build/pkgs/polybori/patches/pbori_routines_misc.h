// -*- c++ -*-
//*****************************************************************************
/** @file pbori_routines_misc.h
 *
 * @author Alexander Dreyer
 * @date 2006-08-23
 *
 * This file includes files, which defines miscellaneous function templates.
 *
 * @par Copyright:
 *   (c) 2006 by The PolyBoRi Team
**/
//*****************************************************************************

// include basic definitions
#include "pbori_defs.h"

// temprarily
#include "CIdxVariable.h"

// temprarily
#include "CacheManager.h"

#include "CDDOperations.h"

BEGIN_NAMESPACE_PBORI

template<class Iterator>
typename Iterator::value_type
index_vector_hash(Iterator start, Iterator finish){

  typedef typename Iterator::value_type value_type;

  value_type vars = 0;
  value_type sum = 0;

  while (start != finish){
    vars++;
    sum += ((*start)+1) * ((*start)+1);
    ++start;
  }
  return sum * vars;
}

/// Function templates for determining the degree of a decision diagram
/// with the help of cache (e. g. CDegreeCache)
template <class DegreeCacher, class NaviType>
typename NaviType::deg_type
dd_cached_degree(const DegreeCacher& cache, NaviType navi) {

  typedef typename NaviType::deg_type deg_type;

  if (navi.isConstant()){ // No need for caching of constant nodes' degrees
    if (navi.terminalValue())
        return 0;
    else
        return -1;
  }

  // Look whether result was cached before
  typename DegreeCacher::node_type result = cache.find(navi);
  if (result.isValid())
    return (*result);

  // Get degree of then branch (contains at least one valid path)...
  deg_type deg = dd_cached_degree(cache, navi.thenBranch()) + 1;

  // ... combine with degree of else branch
  deg = std::max(deg,  dd_cached_degree(cache, navi.elseBranch()) );

  // Write result to cache
  cache.insert(navi, deg);

  return deg;
}

/// Function templates for determining the degree of a decision diagram
/// with the help of cache (e. g. CDegreeCache)
/// Variant *with* given upper bound
/// Assumming that the bound is valid!
template <class DegreeCacher, class NaviType, class SizeType>
typename NaviType::deg_type
dd_cached_degree(const DegreeCacher& cache, NaviType navi, SizeType bound) {

  typedef typename NaviType::deg_type deg_type;

  if (bound < 0) {
    return -1;
  }

  // No need for caching of constant nodes' degrees
  if (bound == 0) {
    while(!navi.isConstant())
      navi.incrementElse();
  }
  if (navi.isConstant())
    return (navi.terminalValue()? 0: -1);

  // Look whether result was cached before
  typename DegreeCacher::node_type result = cache.find(navi, bound);
  if (result.isValid())
    return *result;

  // Get degree of then branch (contains at least one valid path)...
  deg_type deg = dd_cached_degree(cache, navi.thenBranch(), bound - 1);
  if (deg >= 0) ++deg;

  // ... combine with degree of else branch
  if (bound > deg)              // if deg <= bound, we are already finished
    deg = std::max(deg,  dd_cached_degree(cache, navi.elseBranch(), bound) );

  // Write result to cache
  if (deg >= 0) {
    assert(bound >= 0);
    cache.insert(navi, bound, deg);
  }

  return deg;
}

template <class Iterator, class NameGenerator,
          class Separator, class EmptySetType,
          class OStreamType>
void
dd_print_term(Iterator start, Iterator finish, const NameGenerator& get_name,
              const Separator& sep, const EmptySetType& emptyset,
              OStreamType& os){

  if (start != finish){
    os << get_name(*start);
    ++start;
  }
  else
    os << emptyset();

  while (start != finish){
    os << sep() << get_name(*start);
    ++start;
  }
}

template <class TermType, class NameGenerator,
          class Separator, class EmptySetType,
          class OStreamType>
void
dd_print_term(const TermType& term, const NameGenerator& get_name,
              const Separator& sep, const EmptySetType& emptyset,
              OStreamType& os){
  dd_print_term(term.begin(), term.end(), get_name, sep, emptyset, os);
}


template <class Iterator, class NameGenerator,
          class Separator, class InnerSeparator,
          class EmptySetType, class OStreamType>
void
dd_print_terms(Iterator start, Iterator finish, const NameGenerator& get_name,
               const Separator& sep, const InnerSeparator& innersep,
               const EmptySetType& emptyset, OStreamType& os) {

  if (start != finish){
    dd_print_term(*start, get_name, innersep, emptyset, os);
    ++start;
  }

  while (start != finish){
    os << sep();
    dd_print_term(*start, get_name, innersep, emptyset, os);
    ++start;
  }

}


template <bool use_fast,
          class CacheType, class NaviType, class PolyType>
PolyType
dd_multiply(const CacheType& cache_mgr,
            NaviType firstNavi, NaviType secondNavi, PolyType init){

  // Extract subtypes
  typedef typename PolyType::dd_type dd_type;
  typedef typename NaviType::idx_type idx_type;
  typedef NaviType navigator;

  if (firstNavi.isConstant()) {
    if(firstNavi.terminalValue())
      return cache_mgr.generate(secondNavi);
    else
      return cache_mgr.zero();
  }

  if (secondNavi.isConstant()) {
    if(secondNavi.terminalValue())
      return cache_mgr.generate(firstNavi);
    else
      return cache_mgr.zero();
  }
  if (firstNavi == secondNavi)
    return cache_mgr.generate(firstNavi);

  // Look up, whether operation was already used
  navigator cached = cache_mgr.find(firstNavi, secondNavi);
  PolyType result;

  if (cached.isValid()) {       // Cache lookup sucessful
    return cache_mgr.generate(cached);
  }
  else {                        // Cache lookup not sucessful
    // Get top variable's index

    if (*secondNavi < *firstNavi)
      std::swap(firstNavi, secondNavi);

    idx_type index = *firstNavi;

    // Get then- and else-branches wrt. current indexed variable
    navigator as0 = firstNavi.elseBranch();
    navigator as1 = firstNavi.thenBranch();

    navigator bs0;
    navigator bs1;

    if (*secondNavi == index) {
      bs0 = secondNavi.elseBranch();
      bs1 = secondNavi.thenBranch();
    }
    else {
      bs0 = secondNavi;
      bs1 = cache_mgr.zero().navigation();
    }
    PolyType result0 = dd_multiply<use_fast>(cache_mgr, as0, bs0, init);
    PolyType result1;

    // use fast multiplication
    if (use_fast && (*firstNavi == *secondNavi)) {

      PolyType res10 = PolyType(cache_mgr.generate(as1)) +
        PolyType(cache_mgr.generate(as0));
      PolyType res01 = PolyType(cache_mgr.generate(bs0)) +
        PolyType(cache_mgr.generate(bs1));

      result1 = dd_multiply<use_fast>(cache_mgr, res10.navigation(),
                                        res01.navigation(), init) - result0;
    }
    // not using fast multiplication
    else if (as0 == as1) {
      result1 = dd_multiply<use_fast>(cache_mgr, bs0, as1, init);

    }
    else {
      result1 = dd_multiply<use_fast>(cache_mgr, as0, bs1, init);
      if (bs0 != bs1){
        PolyType bs01 = PolyType(cache_mgr.generate(bs0)) +
          PolyType(cache_mgr.generate(bs1));

        result1 +=
          dd_multiply<use_fast>(cache_mgr, bs01.navigation(), as1, init);
      }
    }
    result = dd_type(index, result1.diagram(), result0.diagram());

    // Insert in cache
    cache_mgr.insert(firstNavi, secondNavi, result.navigation());
  } // end of Cache lookup not sucessful

  return result;
}

template <class CacheType, class NaviType, class PolyType>
PolyType
dd_multiply_recursively(const CacheType& cache_mgr,
                        NaviType firstNavi, NaviType secondNavi, PolyType init){

  enum { use_fast =
#ifdef PBORI_FAST_MULTIPLICATION
         true
#else
         false
#endif
  };

  return dd_multiply<use_fast>(cache_mgr, firstNavi, secondNavi, init);
}

template <class CacheType, class NaviType, class PolyType>
PolyType
dd_multiply_recursively_monom(const CacheType& cache_mgr,
                        NaviType monomNavi, NaviType navi, PolyType init){

  // Extract subtypes
  typedef typename PolyType::dd_type dd_type;
  typedef typename NaviType::idx_type idx_type;
  typedef NaviType navigator;

  if (monomNavi.isConstant()) {
    if(monomNavi.terminalValue())
      return cache_mgr.generate(navi);
    else
      return cache_mgr.zero();
  }

  assert(monomNavi.elseBranch().isEmpty());

  if (navi.isConstant()) {
    if(navi.terminalValue())
      return cache_mgr.generate(monomNavi);
    else
      return cache_mgr.zero();
  }
  if (monomNavi == navi)
    return cache_mgr.generate(monomNavi);

  // Look up, whether operation was already used
  navigator cached = cache_mgr.find(monomNavi, navi);

  if (cached.isValid()) {       // Cache lookup sucessful
    return cache_mgr.generate(cached);
  }

  // Cache lookup not sucessful
  // Get top variables' index

  idx_type index = *navi;
  idx_type monomIndex = *monomNavi;

  if (monomIndex < index) {     // Case: index may occure within monom
    init = dd_multiply_recursively_monom(cache_mgr, monomNavi.thenBranch(), navi,
                                   init).diagram().change(monomIndex);
  }
  else if (monomIndex == index) { // Case: monom and poly start with same index

    // Increment navigators
    navigator monomThen = monomNavi.thenBranch();
    navigator naviThen = navi.thenBranch();
    navigator naviElse = navi.elseBranch();

    if (naviThen != naviElse)
      init = (dd_multiply_recursively_monom(cache_mgr, monomThen, naviThen, init)
              + dd_multiply_recursively_monom(cache_mgr, monomThen, naviElse,
                                              init)).diagram().change(index);
  }
  else {                        // Case: var(index) not part of monomial

    init =
      dd_type(index,
              dd_multiply_recursively_monom(cache_mgr, monomNavi, navi.thenBranch(),
                                      init).diagram(),
              dd_multiply_recursively_monom(cache_mgr, monomNavi, navi.elseBranch(),
                                      init).diagram() );
  }

  // Insert in cache
  cache_mgr.insert(monomNavi, navi, init.navigation());

  return init;
}


template <class DDGenerator, class Iterator, class NaviType, class PolyType>
PolyType
dd_multiply_recursively_exp(const DDGenerator& ddgen,
                            Iterator start, Iterator finish,
                            NaviType navi, PolyType init){
  // Extract subtypes
  typedef typename NaviType::idx_type idx_type;
  typedef typename PolyType::dd_type dd_type;
  typedef NaviType navigator;

  if (start == finish)
    return ddgen.generate(navi);

  PolyType result;
  if (navi.isConstant()) {
    if(navi.terminalValue()) {

      std::reverse_iterator<Iterator> rstart(finish), rfinish(start);
      result = ddgen.generate(navi);
      while (rstart != rfinish) {
        result = result.diagram().change(*rstart);
        ++rstart;
      }
    }
    else
      return ddgen.zero();

    return result;
  }

  // Cache lookup not sucessful
  // Get top variables' index

  idx_type index = *navi;
  idx_type monomIndex = *start;

  if (monomIndex < index) {     // Case: index may occure within monom

    Iterator next(start);
    while( (next != finish) && (*next < index) )
      ++next;

    result = dd_multiply_recursively_exp(ddgen, next, finish, navi, init);

    std::reverse_iterator<Iterator> rstart(next), rfinish(start);
    while (rstart != rfinish) {
      result = result.diagram().change(*rstart);
      ++rstart;
    }
  }
  else if (monomIndex == index) { // Case: monom and poly start with same index

    // Increment navigators
    ++start;

    navigator naviThen = navi.thenBranch();
    navigator naviElse = navi.elseBranch();

    if (naviThen != naviElse)
      result =(dd_multiply_recursively_exp(ddgen, start, finish, naviThen, init)
              + dd_multiply_recursively_exp(ddgen, start, finish, naviElse,
                                            init)).diagram().change(index);
  }
  else {                        // Case: var(index) not part of monomial

    result =
      dd_type(index,
              dd_multiply_recursively_exp(ddgen, start, finish,
                                          navi.thenBranch(), init).diagram(),
              dd_multiply_recursively_exp(ddgen, start, finish,
                                          navi.elseBranch(), init).diagram() );
  }

  return result;
}

template<class DegCacheMgr, class NaviType, class SizeType>
bool max_degree_on_then(const DegCacheMgr& deg_mgr, NaviType navi,
                        SizeType degree, valid_tag is_descending) {

  navi.incrementThen();
  return ((dd_cached_degree(deg_mgr, navi, degree - 1) + 1) == degree);
}

template<class DegCacheMgr, class NaviType, class SizeType>
bool max_degree_on_then(const DegCacheMgr& deg_mgr, NaviType navi,
                        SizeType degree, invalid_tag non_descending) {

  navi.incrementElse();
  return (dd_cached_degree(deg_mgr, navi, degree) != degree);
}

template <class NaviType>
NaviType dd_get_constant(NaviType navi) {
  while(!navi.isConstant()) {
    navi.incrementElse();
  }

  return navi;
}

template <class T> class CIndexCacheHandle;
// with degree bound
template <class CacheType, class DegCacheMgr, class NaviType,
          class TermType, class SizeType, class DescendingProperty>
TermType
dd_recursive_degree_lead(const CacheType& cache_mgr, const DegCacheMgr&
                         deg_mgr,
                         NaviType navi, TermType init, SizeType degree,
                         DescendingProperty prop) {


  typedef CIndexCacheHandle<NaviType> deg2node;

  if UNLIKELY(degree < 0)
    return cache_mgr.zero();

  if (navi.isConstant())
    return cache_mgr.generate(navi);

  if (degree == 0)
    return cache_mgr.generate(dd_get_constant(navi));

  // Check cache for previous results
  NaviType cached = cache_mgr.find(navi, deg2node(degree, cache_mgr.manager()));
  if (cached.isValid())
    return cache_mgr.generate(cached);

  // Go to next branch
  if ( max_degree_on_then(deg_mgr, navi, degree, prop) ) {
    NaviType then_branch = navi.thenBranch();
    init = dd_recursive_degree_lead(cache_mgr, deg_mgr, then_branch,
        init, degree - 1, prop);

    if  ((navi.elseBranch().isEmpty()) && (init.navigation()==then_branch))
      init = cache_mgr.generate(navi);
    else
      init = init.change(*navi);

  }
  else {
    init = dd_recursive_degree_lead(cache_mgr, deg_mgr, navi.elseBranch(),
                                    init, degree, prop);
  }
  /// caching deactivated -  wrong results might pollute the cache @todo fix
  NaviType resultNavi(init.navigation());
  cache_mgr.insert(navi, deg2node(degree, cache_mgr.manager()), resultNavi);
//   if (!cache_mgr.find(resultNavi).isValid())
//     deg_mgr.insert(resultNavi, degree, degree);

  return init;
}

template <class CacheType, class DegCacheMgr, class NaviType,
          class TermType, class DescendingProperty>
TermType
dd_recursive_degree_lead(const CacheType& cache_mgr, const DegCacheMgr& deg_mgr,
                         NaviType navi, TermType init, DescendingProperty prop){

//   if (navi.isConstant())
//     return cache_mgr.generate(dd_get_constant(navi));

  return dd_recursive_degree_lead(cache_mgr, deg_mgr, navi, init,
                                  dd_cached_degree(deg_mgr, navi), prop);
}

template <class CacheType, class DegCacheMgr, class NaviType,
          class TermType, class SizeType, class DescendingProperty>
TermType&
dd_recursive_degree_leadexp(const CacheType& cache_mgr,
                            const DegCacheMgr& deg_mgr,
                            NaviType navi, TermType& result,
                            SizeType degree,
                            DescendingProperty prop) {
  typedef CIndexCacheHandle<NaviType> deg2node;
  if UNLIKELY(degree < 0) {
    result.resize(0);
    return result;
  }
  if (navi.isConstant())
    return result;

  if (degree == 0) {
    while (!navi.isConstant())
      navi.incrementElse();
    if (!navi.terminalValue())
      result.resize(0);

    return result;
  }

  // Check cache for previous result
  NaviType cached = cache_mgr.find(navi, deg2node(degree, cache_mgr.manager()));
  if (cached.isValid())
    return result = result.multiplyFirst(cache_mgr.generate(cached));

  // Prepare next branch
  if ( max_degree_on_then(deg_mgr, navi, degree, prop) ) {
    result.push_back(*navi);
    navi.incrementThen();
    --degree;
  }
  else
    navi.incrementElse();

  return
    dd_recursive_degree_leadexp(cache_mgr, deg_mgr, navi, result, degree, prop);
}

template <class CacheType, class DegCacheMgr, class NaviType,
          class TermType, class DescendingProperty>
TermType&
dd_recursive_degree_leadexp(const CacheType& cache_mgr,
                            const DegCacheMgr& deg_mgr,
                            NaviType navi, TermType& result,
                            DescendingProperty prop) {

  if (navi.isConstant())
    return result;

  return dd_recursive_degree_leadexp(cache_mgr, deg_mgr, navi, result,
                                     dd_cached_degree(deg_mgr, navi), prop);
}

// Existential Abstraction. Given a ZDD F, and a set of variables
// S, remove all occurrences of s in S from any subset in F. This can
// be implemented by cofactoring F with respect to s = 1 and s = 0,
// then forming the union of these results.


template <class CacheType, class NaviType, class TermType>
TermType
dd_existential_abstraction(const CacheType& cache_mgr,
                           NaviType varsNavi, NaviType navi, TermType init){

  typedef typename TermType::dd_type dd_type;
  typedef typename TermType::idx_type idx_type;

  if (navi.isConstant())
    return cache_mgr.generate(navi);

  idx_type index(*navi);
  while (!varsNavi.isConstant() && ((*varsNavi) < index))
    varsNavi.incrementThen();

  if (varsNavi.isConstant())
    return cache_mgr.generate(navi);

  // Check cache for previous result
  NaviType cached = cache_mgr.find(varsNavi, navi);
  if (cached.isValid())
    return cache_mgr.generate(cached);

  NaviType thenNavi(navi.thenBranch()), elseNavi(navi.elseBranch());

  TermType thenResult =
    dd_existential_abstraction(cache_mgr, varsNavi, thenNavi, init);

  TermType elseResult =
    dd_existential_abstraction(cache_mgr, varsNavi, elseNavi, init);

  if ((*varsNavi) == index)
    init = thenResult.unite(elseResult);
  else if ( (thenResult.navigation() == thenNavi) &&
            (elseResult.navigation() == elseNavi)  )
    init = cache_mgr.generate(navi);
  else
    init = dd_type(index, thenResult, elseResult);

  // Insert result to cache
  cache_mgr.insert(varsNavi, navi, init.navigation());

  return init;
}



template <class CacheType, class NaviType, class PolyType>
PolyType
dd_divide_recursively(const CacheType& cache_mgr,
                      NaviType navi, NaviType monomNavi, PolyType init){

  // Extract subtypes
  typedef typename NaviType::idx_type idx_type;
  typedef NaviType navigator;
  typedef typename PolyType::dd_type dd_type;

  if (monomNavi.isConstant()) {
    assert(monomNavi.terminalValue() == true);
    return cache_mgr.generate(navi);
  }

  assert(monomNavi.elseBranch().isEmpty());

  if (navi.isConstant())
    return cache_mgr.zero();

  if (monomNavi == navi)
    return cache_mgr.one();

  // Look up, whether operation was already used
  navigator cached = cache_mgr.find(navi, monomNavi);

  if (cached.isValid()) {       // Cache lookup sucessful
    return cache_mgr.generate(cached);
  }

  // Cache lookup not sucessful
  // Get top variables' index

  idx_type index = *navi;
  idx_type monomIndex = *monomNavi;

  if (monomIndex == index) {    // Case: monom and poly start with same index

    // Increment navigators
    navigator monomThen =  monomNavi.thenBranch();
    navigator naviThen = navi.thenBranch();

    init = dd_divide_recursively(cache_mgr, naviThen, monomThen, init);
  }
  else if (monomIndex > index) { // Case: monomIndex may occure within poly

    init =
      dd_type(index,
              dd_divide_recursively(cache_mgr,  navi.thenBranch(), monomNavi,
                                      init).diagram(),
              dd_divide_recursively(cache_mgr, navi.elseBranch(), monomNavi,
                                      init).diagram() );
  }

  // Insert in cache
  cache_mgr.insert(navi, monomNavi,  init.navigation());

  return init;
}



template <class DDGenerator, class Iterator, class NaviType, class PolyType>
PolyType
dd_divide_recursively_exp(const DDGenerator& ddgen,
                          NaviType navi, Iterator start, Iterator finish,
                          PolyType init){

  // Extract subtypes
  typedef typename NaviType::idx_type idx_type;
  typedef typename PolyType::dd_type dd_type;
  typedef NaviType navigator;

  if (start == finish)
    return ddgen.generate(navi);

  if (navi.isConstant())
    return ddgen.zero();


  // Cache lookup not sucessful
  // Get top variables' index

  idx_type index = *navi;
  idx_type monomIndex = *start;

  PolyType result;
  if (monomIndex == index) {    // Case: monom and poly start with same index

    // Increment navigators
    ++start;
    navigator naviThen = navi.thenBranch();

    result = dd_divide_recursively_exp(ddgen, naviThen, start, finish, init);
  }
  else if (monomIndex > index) { // Case: monomIndex may occure within poly

    result =
      dd_type(index,
              dd_divide_recursively_exp(ddgen, navi.thenBranch(), start, finish,
                                        init).diagram(),
              dd_divide_recursively_exp(ddgen, navi.elseBranch(), start, finish,
                                        init).diagram() );
  }
  else
    result = ddgen.zero();

  return result;
}

/// Function templates for determining the used variables of a decision diagram
/// with the help of cache
template <class CacheType, class NaviType, class MonomType>
MonomType
cached_used_vars(const CacheType& cache, NaviType navi, MonomType init) {

  if (navi.isConstant()) // No need for caching of constant nodes' degrees
    return init;

  // Look whether result was cached before
  NaviType cached_result = cache.find(navi);

  typedef typename MonomType::poly_type poly_type;
  if (cached_result.isValid())
    return CDDOperations<typename MonomType::dd_type,
    MonomType>().getMonomial(cache.generate(cached_result));

  MonomType result = cached_used_vars(cache, navi.thenBranch(), init);
  result *= cached_used_vars(cache, navi.elseBranch(), init);

  result = result.change(*navi);

  // Write result to cache
  cache.insert(navi, result.diagram().navigation());

  return result;
}

template <class NaviType, class Iterator>
bool
dd_owns(NaviType navi, Iterator start, Iterator finish) {

  if (start == finish) {
    while(!navi.isConstant())
      navi.incrementElse();
    return navi.terminalValue();
  }

  while(!navi.isConstant() && (*start > *navi) )
    navi.incrementElse();

  if (navi.isConstant() || (*start != *navi))
    return false;

  return dd_owns(navi.thenBranch(), ++start, finish);
}

/// Test whether a set of monomials contains all divisors of degree - 1 of a
/// given monomial
template <class NaviType, class MonomIterator>
bool
dd_contains_divs_of_dec_deg(NaviType navi,
                            MonomIterator start, MonomIterator finish) {

  // Managing trivial cases

  if (start == finish)          // divisors of monomial 1 is the empty set,
                                // which is always contained
    return true;

  if (navi.isConstant()) {      // set is empty or owns 1 only
    if (navi.terminalValue())   // set == {1}
      return (++start == finish);  // whether monom is of degree one
    else                        // empty set contains no divisors
      return false;

  }

  // Actual computations

  if (*navi < *start)
    return dd_contains_divs_of_dec_deg(navi.elseBranch(), start, finish);


  if (*navi > *start) {
    if (++start == finish)      // if monom is of degree one
      return owns_one(navi);
    else
      return false;
  }

  // (*navi == *start) here
  ++start;
  return dd_owns(navi.elseBranch(), start, finish) &&
    dd_contains_divs_of_dec_deg(navi.thenBranch(), start, finish);
}



// determine the part of a polynomials of a given degree
template <class CacheType, class NaviType, class DegType, class SetType>
SetType
dd_graded_part(const CacheType& cache, NaviType navi, DegType deg,
               SetType init) {


  if (deg == 0) {
    while(!navi.isConstant())
      navi.incrementElse();
    return cache.generate(navi);
  }

  if(navi.isConstant())
    return cache.zero();

  // Look whether result was cached before
  NaviType cached = cache.find(navi, deg);

  if (cached.isValid())
    return cache.generate(cached);

  SetType result =
    SetType(*navi,
            dd_graded_part(cache, navi.thenBranch(), deg - 1, init),
            dd_graded_part(cache, navi.elseBranch(), deg, init)
            );

  // store result for later reuse
  cache.insert(navi, deg, result.navigation());

  return result;
}


/// Function templates extracting the terms of a given decision diagram contain
/// which contains only indices from first lexicographical path in
/// Note: Replacement for dd_intersect_some_index
template <class CacheManager, class NaviType, class SetType>
SetType
dd_first_divisors_of(CacheManager cache_mgr, NaviType navi,
                     NaviType rhsNavi, SetType init ) {

  typedef typename SetType::dd_type dd_type;
  while( (!navi.isConstant()) && (*rhsNavi != *navi) ) {

    if ( (*rhsNavi < *navi) && (!rhsNavi.isConstant()) )
      rhsNavi.incrementThen();
    else
      navi.incrementElse();
  }

  if (navi.isConstant())        // At end of path
    return cache_mgr.generate(navi);

  // Look up, whether operation was already used
  NaviType result = cache_mgr.find(navi, rhsNavi);

  if (result.isValid())       // Cache lookup sucessful
    return  cache_mgr.generate(result);

  assert(*rhsNavi == *navi);

  // Compute new result
  init = dd_type(*rhsNavi,
                 dd_first_divisors_of(cache_mgr, navi.thenBranch(), rhsNavi,
                                      init).diagram(),
                 dd_first_divisors_of(cache_mgr, navi.elseBranch(), rhsNavi,
                                      init).diagram() );
  // Insert result to cache
  cache_mgr.insert(navi, rhsNavi, init.navigation());

  return init;
}

template <class CacheType, class NaviType, class SetType>
SetType
dd_first_multiples_of(const CacheType& cache_mgr,
                      NaviType navi, NaviType rhsNavi, SetType init){

  typedef typename SetType::dd_type dd_type;

  if(rhsNavi.isConstant()) {
    assert(rhsNavi.terminalValue() == true);
    return cache_mgr.generate(navi);
  }

  if (navi.isConstant() || (*navi > *rhsNavi))
    return cache_mgr.zero();

  if (*navi == *rhsNavi)
    return dd_first_multiples_of(cache_mgr, navi.thenBranch(),
                                 rhsNavi.thenBranch(), init).change(*navi);

  // Look up old result - if any
  NaviType result = cache_mgr.find(navi, rhsNavi);

  if (result.isValid())
    return cache_mgr.generate(result);

  // Compute new result
  init = dd_type(*navi,
                  dd_first_multiples_of(cache_mgr, navi.thenBranch(),
                                        rhsNavi, init).diagram(),
                  dd_first_multiples_of(cache_mgr, navi.elseBranch(),
                                        rhsNavi, init).diagram() );

  // Insert new result in cache
  cache_mgr.insert(navi, rhsNavi, init.navigation());

  return init;
}


/// Internal variant: Maps a polynomial from one ring to another, using a given
/// map old_index -> new_polynomial
template <class PolyType, class RingType, class MapType, class NaviType>
PolyType
substitute_variables__(const RingType& ring, const MapType& idx2poly, NaviType navi) {

  if (navi.isConstant())
    return ring.constant(navi.terminalValue());

  typename MapType::const_reference var(idx2poly[*navi]);
  //  assert(var.ring() == ring);
  return (idx2poly[*navi] *
          substitute_variables__<PolyType>(ring, idx2poly, navi.thenBranch())) +
    substitute_variables__<PolyType>(ring, idx2poly, navi.elseBranch());
}

/// Maps a polynomial from one ring to another, using a given map
/// old_index -> new_polynomial
template <class RingType, class MapType, class PolyType>
PolyType
substitute_variables(const RingType& ring,
                     const MapType& idx2poly, const PolyType& poly) {

  return substitute_variables__<PolyType>(ring, idx2poly, poly.navigation());
}


END_NAMESPACE_PBORI
