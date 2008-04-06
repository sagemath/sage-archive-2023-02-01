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


    void pbenv_changeOrdering "BooleEnv::changeOrdering"(ordercodes c)
    ordercodes  pbenv_getOrderCode "BooleEnv::getOrderCode"()
    bint pbenv_isDegreeOrder "BooleEnv::ordering().isDegreeOrder"()


    ctypedef struct PBRing "BoolePolyRing":
        int (* nVariables)()
        PBDD (* variable)(int n)
        void (*activate)()

    ctypedef struct PBMonomIter "BooleMonomial::const_iterator":
        int (* value "operator*")()
        int (* next "operator++")()
        int (* hash)()

    void PBMonomIter_destruct "Destruct<BooleMonomial::const_iterator>" \
            (PBMonomIter *mem)

    ctypedef struct PBVar "BooleVariable":
        int (* index)()
        bint (* is_equal "operator==")(PBVar right)

    PBVar* PBVar_construct_int "Construct_p<BooleVariable, int>" \
        (void *mem, int ind)

    PBVar* PBVar_construct_pbvar "Construct_p<BooleVariable, BooleVariable>" \
        (void *mem, PBVar v)

    ctypedef struct PBMonomVarIter "BooleMonomial::variable_iterator":
        PBVar (* value "operator*")()
        int (* next "operator++")()
        bint (* equal "operator==")(PBMonomVarIter rhs)

    void PBMonomVarIter_destruct "Destruct<BooleMonomial::variable_iterator>" \
            (PBMonomVarIter *mem)

    ctypedef struct PBSet "BooleSet"
    ctypedef struct PBMonom "BooleMonomial":
        bint (* reducibleBy)(PBMonom rhs)
        int (* deg)()
        int (* hash)()
        int (* stableHash)()
        int (* firstIndex)()
        comparecodes (* compare)(PBMonom rhs)
        PBSet (* set)()
        PBSet (* divisors)()
        PBSet (* multiples)(PBMonom rhs)
        PBMonomIter (* begin)()
        PBMonomIter (* end)()
        PBMonomVarIter (*variableBegin)()
        PBMonomVarIter (*variableEnd)()
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

    PBSetIter* PBSetIter_construct \
            "Construct<BooleSet::const_iterator>"(PBSetIter *mem)
    void PBSetIter_destruct "Destruct<BooleSet::const_iterator>"(PBSetIter *mem)

    ctypedef struct PBSet "BooleSet":
        bint (* emptiness)()
        bint (* owns)(PBMonom val)
        int (* nNodes)()
        int (* nSupport)()
        int (* length)()
        int (* hash)()
        int (* stableHash)()
        PBNavigator (* navigation)()
        PBSet (* cartesianProduct)(PBSet rhs)
        PBSet (* diff)(PBSet rhs)
        PBSet (* divide)(PBMonom rhs)
        PBSet (* change)(int idx)
        PBSet (* subset0)(int i)
        PBSet (* subset1)(int i)
        PBSet (* unite)(PBSet rhs)
        PBMonom (* usedVariables)()
        PBSetIter (* begin)()
        PBSetIter (* end)()

    PBSet pb_include_divisors "include_divisors" (PBSet p)
    PBSet pb_minimal_elements "minimal_elements" (PBSet p)

    object PBSet_to_str "_to_PyString<BooleSet>"(PBSet *p)

    # non-allocating versions
    PBSet* PBSet_construct "Construct<BooleSet>"(void* mem)
    PBSet* PBSet_construct_pbset \
            "Construct_p<BooleSet, BooleSet>" (void* mem, PBSet d)
    PBSet* PBSet_construct_pbnav \
            "Construct_pp<BooleSet, CCuddNavigator, BooleRing>" \
            (void* mem, PBNavigator d, PBRing r)
    PBSet* PBSet_construct_indsetset \
            "Construct_pppp<BooleSet, int, CCuddNavigator, CCuddNavigator, BooleRing>" \
            (void* mem, int ind, PBNavigator a, PBNavigator b, PBRing r)
    void PBSet_destruct "Destruct<BooleSet>"(PBSet *mem)


    ctypedef struct PBPolyIter "BoolePolynomial::ordered_iterator":
        PBMonom (* value "operator*")()
        int (* next "operator++")()
        bint (* equal)(PBPolyIter rhs)

    void PBPolyIter_destruct "Destruct<BoolePolynomial::ordered_iterator>"(PBPolyIter *mem)

    ctypedef struct PBPoly "BoolePolynomial":
        int (* deg)()
        int (* lmDeg)()
        int (* lexLmDeg)()
        int (* totalDeg)()
        int (* length)()
        int (* eliminationLength)()
        int (* nUsedVariables)()
        int (* nNodes)()
        int (* hash)()
        int (* stableHash)()
        bint (* isZero)()
        bint (* isOne)()
        bint (* isConstant)()
        bint (* isSingleton)()
        bint (* hasConstantPart)()
        bint (* reducibleBy)(PBPoly rhs)
        PBMonom (* lead)()
        PBMonom (* lexLead)()
        PBMonom (* firstTerm)()
        PBMonom (* usedVariables)()
        PBPoly (* gradedPart)(int d)
        PBDD (* diagram)()
        PBSet (* set)()
        PBSet (* lmDivisors)()
        PBNavigator (* navigation)()
        PBPolyIter (* orderedBegin)()
        PBPolyIter (* orderedEnd)()
        void (* iadd "operator+=")(PBPoly right)
        void (* iadd_PBMonom "operator+=")(PBMonom right)
        void (* imul "operator*=")(PBPoly right)
        void (* imul_monom "operator*=")(PBMonom right)
        bint (* is_equal "operator==")(PBPoly right)

    PBSet pb_zeroes "zeroes" (PBPoly p, PBSet s)
    PBPoly pb_spoly "spoly" (PBPoly p, PBPoly r)

    PBPoly pb_map_every_x_to_x_plus_one "map_every_x_to_x_plus_one" (PBPoly)

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
        bint (* variableHasValue)(int idx)
        int (* addGenerator)(PBPoly, bint is_impl)
        void (* addGeneratorDelayed)(PBPoly)
        void (* addAsYouWish)(PBPoly)
        void (* symmGB_F2)()
        void (* llReduceAll)()
        void (* cleanTopByChainCriterion "pairs.cleanTopByChainCriterion")()
        int (* nGenerators "generators.size")()
        int (* npairs "pairs.queue.size")()
        int (* suggestPluginVariable)()
        PBPoly (* nextSpoly)()
        PBPoly (* nf)(PBPoly p)
        PBPolyVector (* allGenerators)()
        PBPolyVector (* minimalize)()
        PBPolyVector (* minimalizeAndTailReduce)()
        PBPolyVector (* faugereStepDense)(PBPolyVector v)

    int (* pairs_top_sugar)(GBStrategy strat)
    PBPolyVector (* someNextDegreeSpolys)(GBStrategy strat, int n)
    PBPolyVector (* nextDegreeSpolys)(GBStrategy strat)
    PBPolyVector (* small_next_degree_spolys)(GBStrategy strat, double f, int n)
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

    PBRing get_current_ring "BooleEnv::ring"()

    PBPoly pb_add_up_polynomials "add_up_polynomials" \
            (PBPolyVector v)

    PBPoly pb_nf3 "nf3" (GBStrategy strat, PBPoly p, PBMonom m)

    PBPoly pb_red_tail "red_tail" (GBStrategy strat, PBPoly p)

    PBPoly pb_interpolate "interpolate" (PBSet z, PBSet o)

    PBPoly pb_interpolate_smallest_lex "interpolate_smallest_lex" \
            (PBSet z, PBSet o)

    PBSet pb_contained_variables_cudd_style "contained_variables_cudd_style" \
            (PBSet m)

    PBSet pb_mod_var_set "mod_var_set" (PBSet a, PBSet v)

    PBPoly pb_mult_fast_sim "mult_fast_sim" (PBPolyVector v)

    int pb_select1 "select1" (GBStrategy s, PBMonom m)

    PBSet pb_recursively_insert "recursively_insert"(PBNavigator p,
                                                int idx, PBSet m)
    PBPoly pb_ll_red_nf_noredsb "ll_red_nf_noredsb"(PBPoly p,
                                                PBSet reductors)
    PBPoly pb_ll_red_nf "ll_red_nf"(PBPoly p, PBSet reductors)

    PBSet pb_mod_mon_set "mod_mon_set"(PBSet as, PBSet vs)

    PBPolyVector pb_parallel_reduce "parallel_reduce" \
        (PBPolyVector inp, GBStrategy strat, int average_steps, double delay_f)

    void pb_append_block "BooleEnv::appendBlock" (int ind)

    void pb_set_variable_name "BooleEnv::setVariableName" \
        (int idx, char *varname)

    #M4RI initialization
    void buildAllCodes()
    void destroyAllCodes()
    void setupPackingMasks()
