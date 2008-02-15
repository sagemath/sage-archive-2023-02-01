cdef extern from "pb_wrap.h":
    cdef enum ordercodes "COrderEnums::ordercodes":
        lp            "CTypes::lp"
        dlex          "CTypes::dlex"
        dp_asc        "CTypes::dp_asc"
        block_dlex    "CTypes::block_dlex"
        block_dp_asc  "CTypes::block_dp_asc"

    cdef enum comparecodes "CCompareEnums::comparecodes":
        less_than               "CTypes::less_than"
        equality                "CTypes::equality"
        greater_than            "CTypes::greater_than"
        less_or_equal_max       "CTypes::less_or_equal_max"
        greater_or_equal_min    "CTypes::greater_or_equal_min"

    ctypedef struct PBNavigator "struct CCuddNavigator":
        PBNavigator (* thenBranch) ()
        PBNavigator (* elseBranch) ()
        int (* value "operator*")()
        bint (* isConstant)()
        bint (* isTerminated)()

    # non-allocating versions
    PBNavigator* PBNavigator_construct \
            "Construct_p<CCuddNavigator, CCuddNavigator>" \
            (void* mem, PBNavigator N)
    void PBNavigator_destruct "Destruct<CCuddNavigator>"(PBNavigator *mem)

    ctypedef struct PBDD "struct CDDInterface<CTypes::dd_base>":
        bint (* emptiness)()
        PBNavigator (* navigation)()
        PBDD (* subset0)(int idx)
        PBDD (* subset1)(int idx)
        PBDD (* unite)(PBDD rhs)

    # non-allocating versions
    void PBDD_destruct "Destruct<CDDInterface<CTypes::dd_base> >"(PBDD *mem)

    ctypedef struct PBRing "BoolePolyRing":
        bint (* isDegreeOrder)()
        int (* nVariables)()
        PBDD (* variable)(int n)
        void (*setRingVariableName)(int idx, char *varname)
        void (*activate)()
        ordercodes (* getOrderCode)()

    ctypedef struct PBMonomIter "BooleMonomial::const_iterator":
        int (* value "operator*")()
        int (* next "operator++")()
        int (* hash)()

    void PBMonomIter_destruct "Destruct<BooleMonomial::const_iterator>" (PBMonomIter *mem)

    ctypedef struct PBVar "BooleVariable":
        int (* index)()
        bint (* is_equal "operator==")(PBVar right)

    PBVar* PBVar_construct_int "Construct_p<BooleVariable, int>" \
        (void *mem, int ind)

    PBVar* PBVar_construct_pbvar "Construct_p<BooleVariable, BooleVariable>" \
        (void *mem, PBVar v)

    ctypedef struct PBMonom "BooleMonomial":
        int (* deg)()
        int (* hash)()
        comparecodes (* compare)(PBMonom rhs)
        PBMonomIter (* begin)()
        PBMonomIter (* end)()
        void (* imul "operator*=")(PBMonom right)

    object PBMonom_to_str "_to_PyString<BooleMonomial>"(PBMonom *p)

    # non-allocating versions
    PBMonom* PBMonom_construct "Construct<BooleMonomial>"(void *mem)
    PBMonom* PBMonom_construct_pbmonom \
            "Construct_p<BooleMonomial, BooleMonomial>" (void *mem, PBMonom m)
    PBMonom* PBMonom_construct_pbvar \
            "Construct_p<BooleMonomial, BooleVariable>" (void *mem, PBVar m)
    PBMonom* PBMonom_construct_dd \
            "Construct_p<BoolePolynomial, BooleMonomial::dd_type>" (void *mem, PBDD m)
    void PBMonom_destruct "Destruct<BooleMonomial>"(PBMonom *mem)

    ctypedef struct PBSetIter "BooleSet::const_iterator":
        PBMonom (* value "operator*")()
        int (* next "operator++")()
        bint (* equal)(PBSetIter rhs)

    void PBSetIter_destruct "Destruct<BooleSet::const_iterator>"(PBSetIter *mem)

    ctypedef struct PBSet "BooleSet":
        bint (* emptiness)()
        bint (* owns)(PBMonom val)
        PBNavigator (* navigation)()
        PBSet (* cartesianProduct)(PBSet rhs)
        PBSet (* diff)(PBSet diff)
        PBSet (* change)(int idx)
        PBSet (* unite)(PBSet rhs)
        PBMonom (* usedVariables)()
        PBSetIter (* begin)()
        PBSetIter (* end)()

    object PBSet_to_str "_to_PyString<BooleSet>"(PBSet *p)

    # non-allocating versions
    PBSet* PBSet_construct "Construct<BooleSet>"(void* mem)
    PBSet* PBSet_construct_pbset \
            "Construct_p<BooleSet, BooleSet>" (void* mem, PBSet d)
    PBSet* PBSet_construct_dd \
            "Construct_p<BooleSet, BooleSet::dd_type>" (void* mem, PBDD d)
    PBSet* PBSet_construct_pbnav \
            "Construct_p<BooleSet, CCuddNavigator>" (void* mem, PBNavigator d)
    PBSet* PBSet_construct_indsetset \
            "Construct_ppp<BooleSet, int, CCuddNavigator, CCuddNavigator>" \
            (void* mem, int ind, PBNavigator a, PBNavigator b)
    void PBSet_destruct "Destruct<BooleSet>"(PBSet *mem)


    ctypedef struct PBPolyIter "BoolePolynomial::ordered_iterator":
        PBMonom (* value "operator*")()
        int (* next "operator++")()
        bint (* equal)(PBPolyIter rhs)

    void PBPolyIter_destruct "Destruct<BoolePolynomial::ordered_iterator>"(PBPolyIter *mem)

    ctypedef struct PBPoly "BoolePolynomial":
        int (* deg)()
        int (* lmDeg)()
        int (* length)()
        int (* eliminationLength)()
        int (* nUsedVariables)()
        bint (* isZero)()
        bint (* isOne)()
        bint (* isConstant)()
        bint (* isSingleton)()
        PBMonom (* lead)()
        PBMonom (* lexLead)()
        PBMonom (* usedVariables)()
        PBDD (* diagram)()
        PBSet (* set)()
        PBNavigator (* navigation)()
        PBPolyIter (* orderedBegin)()
        PBPolyIter (* orderedEnd)()
        void (* iadd "operator+=")(PBPoly right)
        void (* iadd_PBMonom "operator+=")(PBMonom right)
        void (* imul "operator*=")(PBPoly right)
        void (* imul_monom "operator*=")(PBMonom right)
        bint (* is_equal "operator==")(PBPoly right)

    PBPoly map_every_x_to_x_plus_one(PBPoly)

    # non-allocating versions
    PBRing* PBRing_construct \
            "Construct_pp<BoolePolyRing, BoolePolyRing::size_type, BoolePolyRing::ordercode_type> " (void *mem, int nvars, ordercodes order)
    void PBRing_destruct "Destruct<BoolePolyRing>"(PBRing *mem)

    # non-allocating versions
    PBPoly* PBPoly_construct "Construct<BoolePolynomial>"(void *mem)
    PBPoly* PBPoly_construct_dd \
            "Construct_p<BoolePolynomial, BoolePolyRing::dd_type>" \
            (void *mem, PBDD d)
    PBPoly* PBPoly_construct_pbset \
            "Construct_p<BoolePolynomial, BooleSet>" (void *mem, PBSet d)
    PBPoly* PBPoly_construct_pbpoly \
            "Construct_p<BoolePolynomial, BoolePolynomial>" \
            (void *mem, PBPoly d)
    PBPoly* PBPoly_construct_pbmonom \
            "Construct_p<BoolePolynomial, BooleMonomial>" (void *mem, PBMonom d)
    PBPoly* PBPoly_construct_int \
            "Construct_p<BoolePolynomial, int>" (void *mem, int d)
    void PBPoly_destruct "Destruct<BoolePolynomial>"(PBPoly *mem)

    object PBPoly_to_str "_to_PyString<BoolePolynomial>"(PBPoly *p)

    ctypedef struct PBPolyVectorIter \
            "std::vector<BoolePolynomial>::iterator ":
        PBPoly (* value "operator*")()
        int (* next "operator++")()

    void PBPolyVectorIter_destruct "Destruct<std::vector<BoolePolynomial>::iterator>"(PBPolyVectorIter *mem)

    bint PBPolyVectorIter_equal "operator=="(PBPolyVectorIter lhs, \
            PBPolyVectorIter rhs)

    ctypedef struct PBPolyVector "std::vector<BoolePolynomial>":
        int (* size)()
        PBPoly (* get "operator[]")(int)
        PBPolyVectorIter (* begin)()
        PBPolyVectorIter (* end)()
        void (* push_back)(PBPoly val)

    PBPolyVector* PBPolyVector_construct \
            "Construct< std::vector<BoolePolynomial> >"(void *mem)
    void PBPolyVector_destruct "Destruct< std::vector<BoolePolynomial> >"\
            (PBPolyVector *mem)

    ctypedef struct GBStrategy "GroebnerStrategy":
        bint reduceByTailReduced
        bint enabledLog
        unsigned int reductionSteps
        int normalForms
        int currentDegree
        int chainCriterions
        int variableChainCriterions
        int easyProductCriterions
        int extendedProductCriterions
        int averageLength
        bint optRedTail
        bint optLazy
        bint optLL
        bint optDelayNonMinimals
        bint optBrutalReductions
        bint optExchange
        bint optAllowRecursion
        bint optRedTailDegGrowth
        bint optStepBounded
        bint optLinearAlgebraInLastBlock
        bint optRedTailInLastBlock
        PBSet monomials
        PBSet llReductor
        PBSet minimalLeadingTerms
        PBSet leadingTerms
        bint (* containsOne)()
        int (* addGenerator)(PBPoly, bint is_impl)
        void (* addGeneratorDelayed)(PBPoly)
        void (* addAsYouWish)(PBPoly)
        void (* symmGB_F2)()
        void (* cleanTopByChainCriterion "pairs.cleanTopByChainCriterion")()
        int (* nGenerators "generators.size")()
        int (* npairs "pairs.queue.size")()
        PBPolyVector (* minimalize)()
        PBPolyVector (* minimalizeAndTailReduce)()
        PBPolyVector (* faugereStepDense)(PBPolyVector v)

    int (* pairs_top_sugar)(GBStrategy strat)
    PBPolyVector (* someNextDegreeSpolys)(GBStrategy strat, int n)
    void (* implications)(GBStrategy strat, int i)

    PBPoly GB_get_ith_gen "get_ith_gen" (GBStrategy strat, int i)

    # non-allocating versions
    GBStrategy* GBStrategy_construct "Construct<GroebnerStrategy>"(void *mem)
    GBStrategy* GBStrategy_construct_gbstrategy \
            "Construct_p<GroebnerStrategy, GroebnerStrategy>" \
            (void *mem, GBStrategy strat)
    PBPoly* PBPoly_construct_dd \
            "Construct_p<BoolePolynomial, BoolePolyRing::dd_type>" \
            (void *mem, PBDD d)
    void GBStrategy_destruct "Destruct<GroebnerStrategy>"(GBStrategy *mem)

    PBRing get_current_ring "BoolePolyRing::ring"()

    PBSet pb_recursively_insert "recursively_insert"(PBNavigator p,
                                                int idx, PBNavigator m)
    PBPoly pb_ll_red_nf_noredsb "ll_red_nf_noredsb"(PBPoly p,
                                                PBSet reductors)
    PBPoly pb_ll_red_nf "ll_red_nf"(PBPoly p, PBSet reductors)

    PBSet pb_mod_mon_set "mod_mon_set"(PBSet as, PBSet vs)

    void pb_change_ordering "BoolePolyRing::changeOrdering" (ordercodes c)

    PBPolyVector pb_parallel_reduce "parallel_reduce" \
        (PBPolyVector inp, GBStrategy strat, int average_steps, double delay_f)

    void pb_set_variable_name "BoolePolyRing::setRingVariableName" \
        (int idx, char* varname)

    void pb_append_ring_block "BoolePolyRing::appendRingBlock" (int ind)

    #M4RI initialization
    void buildAllCodes()
    void destroyAllCodes()
    void setupPackingMasks()
