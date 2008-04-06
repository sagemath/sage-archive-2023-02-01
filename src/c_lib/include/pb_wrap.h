#include "polybori.h"
#include "groebner_alg.h"
#include "nf.h"
#include "interpolate.h"
#include "ccobject.h"

// M4RI
#define PACKED 1
#include "M4RI/packedmatrix.h"
#include "M4RI/grayflex.h"

#include <sstream>
#include <vector>

USING_NAMESPACE_PBORI
USING_NAMESPACE_PBORIGB

BoolePolynomial get_ith_gen(const GroebnerStrategy& strat, int i){
    return strat.generators[i].p;
}

static int pairs_top_sugar(const GroebnerStrategy& strat){
    if (strat.pairs.pairSetEmpty())
        return -1;
    else
        return (strat.pairs.queue.top().sugar);
}

static std::vector<BoolePolynomial> someNextDegreeSpolys(GroebnerStrategy& strat, int n){
    std::vector<BoolePolynomial> res;
    assert(!(strat.pairs.pairSetEmpty()));
    strat.pairs.cleanTopByChainCriterion();
    deg_type deg=strat.pairs.queue.top().sugar;

    while((!(strat.pairs.pairSetEmpty())) && \
                (strat.pairs.queue.top().sugar<=deg) && (res.size()<n)){
        assert(strat.pairs.queue.top().sugar==deg);
        res.push_back(strat.nextSpoly());
        strat.pairs.cleanTopByChainCriterion();
    }
    return res;
}

static std::vector<Polynomial> nextDegreeSpolys(GroebnerStrategy& strat){
    std::vector<Polynomial> res;
    assert(!(strat.pairs.pairSetEmpty()));
    strat.pairs.cleanTopByChainCriterion();
    deg_type deg=strat.pairs.queue.top().sugar;

    while((!(strat.pairs.pairSetEmpty())) &&
            (strat.pairs.queue.top().sugar<=deg)){

        assert(strat.pairs.queue.top().sugar==deg);
        res.push_back(strat.nextSpoly());
        strat.pairs.cleanTopByChainCriterion();
    }
    return res;
}

static std::vector<Polynomial> small_next_degree_spolys(GroebnerStrategy& strat,
        double f, int n){
    std::vector<Polynomial> res;
    assert(!(strat.pairs.pairSetEmpty()));
    strat.pairs.cleanTopByChainCriterion();
    deg_type deg=strat.pairs.queue.top().sugar;
    wlen_type wlen=strat.pairs.queue.top().wlen;
    while((!(strat.pairs.pairSetEmpty())) &&
            (strat.pairs.queue.top().sugar<=deg) &&
            (strat.pairs.queue.top().wlen<=wlen*f+2)&& (res.size()<n)){

        assert(strat.pairs.queue.top().sugar==deg);
        res.push_back(strat.nextSpoly());
        strat.pairs.cleanTopByChainCriterion();
    }
    return res;
}

static void implications(GroebnerStrategy& strat, int i){
    strat.addNonTrivialImplicationsDelayed(strat.generators[i]);
}

