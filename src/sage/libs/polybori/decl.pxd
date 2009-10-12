cdef extern from "pb_wrap.h":
    ctypedef struct std_string "std::string":
        char *(* c_str)()

    std_string new_stdstring "std::string"(char *str)

    cdef enum ordercodes "COrderEnums::ordercodes":
        pblp            "CTypes::lp"
        pbdlex          "CTypes::dlex"
        pbdp_asc        "CTypes::dp_asc"
        pbblock_dlex    "CTypes::block_dlex"
        pbblock_dp_asc  "CTypes::block_dp_asc"

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
        bint (* is_equal "operator==")(PBNavigator right)
        int (* hash)()

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

    ctypedef struct PBBlockIter "CDynamicOrderBase::block_iterator":
        int (* value "operator*")()
        PBBlockIter (* next "operator++")()
        int (* hash)()

    bint PBBlockIter_equals "operator=="(PBBlockIter lhs, PBBlockIter rhs)

    ctypedef struct PBOrdering "CDynamicOrderBase":
        int (*getOrderCode)()
        int (*getBaseOrderCode)()
        PBBlockIter (*blockBegin)()
        PBBlockIter (*blockEnd)()


    ctypedef struct PBRing "BoolePolyRing":
        int (* nVariables)()
        PBDD (* variable)(int n)
        void (*activate)()
        PBOrdering (*ordering)()
        int (*hash)()
        void (*setVariableName) (int idx, std_string varname)
        std_string (*getVariableName)(int idx)
        PBRing (*clone)()

    void pbenv_changeOrdering "BooleEnv::changeOrdering"(int c)
    int  pbenv_getOrderCode "BooleEnv::getOrderCode"()
    bint pbenv_isDegreeOrder "BooleEnv::ordering().isDegreeOrder"()
    void pbenv_set "BooleEnv::set"(PBRing ring)
    PBRing pbenv_ring "BooleEnv::ring"()

    ctypedef struct PBMonomIter "BooleMonomial::const_iterator":
        int (* value "operator*")()
        int (* next "operator++")()
        bint (* equal "equal")(PBMonomIter rhs)
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
        bint (* equal "equal")(PBMonomVarIter rhs)

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
        void (* idiv "operator/=")(PBMonom right)
        PBNavigator (* navigation)()
        PBMonom (* GCD)(PBMonom rhs)

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
        PBSet (* intersect)(PBSet rhs)
        PBMonom (* usedVariables)()
        PBSet (* multiplesOf)(PBMonom rhs)
        double (* sizeDouble)()
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
        int (* leadDeg)()
        int (* lexLeadDeg)()
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
        PBSet (* leadDivisors)()
        PBNavigator (* navigation)()
        PBPolyIter (* orderedBegin)()
        PBPolyIter (* orderedEnd)()
        void (* iadd "operator+=")(PBPoly right)
        void (* iadd_PBMonom "operator+=")(PBMonom right)
        void (* imul "operator*=")(PBPoly right)
        void (* imul_monom "operator*=")(PBMonom right)
        bint (* is_equal "operator==")(PBPoly right)

    PBSet pb_zeros "zeros" (PBPoly p, PBSet s)
    PBPoly pb_spoly "spoly" (PBPoly p, PBPoly r)

    PBPoly pb_map_every_x_to_x_plus_one "map_every_x_to_x_plus_one" (PBPoly)

    # non-allocating versions
    PBRing* PBRing_construct \
            "Construct_ppp<BoolePolyRing, BoolePolyRing::size_type, BoolePolyRing::ordercode_type, bool> " (void *mem, int nvars, int order, bint make_active)

    PBRing* PBRing_construct_pbring \
            "Construct_p<BoolePolyRing, BoolePolyRing> " (void *mem, PBRing other)

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

    bint PBPolyVectorIter_equal "operator=="(PBPolyVectorIter lhs, PBPolyVectorIter rhs)

    void PBPolyVectorIter_destruct "Destruct<std::vector<BoolePolynomial>::iterator>"(PBPolyVectorIter *mem)

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

    ctypedef struct PBPolyEntry "PolyEntry":
        PBPoly p

    ctypedef struct PBRedStrategy "ReductionStrategy":
        PBSet leadingTerms
        PBSet minimalLeadingTerms
        PBSet leadingTerms11
        PBSet leadingTerms00
        PBSet llReductor
        PBSet monomials
        PBSet monomials_plus_one
        bint optBrutalReductions
        bint optLL
        PBPoly (* nf)(PBPoly p)
        bint optRedTailDegGrowth
        bint optRedTail
        void (* setupSetsForLastElement)()
        bint (* canRewrite)(PBPoly p)
        void (* addGenerator)(PBPoly p)
        int select1(PBMonom p)

        int (* select_short)(PBPoly p)
        PBPoly (* headNormalForm)(PBPoly p)
        PBPoly (* reducedNormalForm)(PBPoly p)

        int (* size)()
        PBPolyEntry (* get "operator[]")(int)

    ctypedef struct PBFglmStrategy "FGLMStrategy":
        PBPolyVector (* main)()

    ctypedef struct PBGBStrategy "GroebnerStrategy":
        bint reduceByTailReduced
        PBRedStrategy generators

        bint enabledLog
        unsigned int reductionSteps
        int normalForms
        int currentDegree
        int chainCriterions
        int variableChainCriterions
        int easyProductCriterions
        int extendedProductCriterions
        int averageLength

        bint optHFE
        bint optLazy

        bint optDelayNonMinimals

        bint optExchange
        bint optAllowRecursion

        bint optStepBounded
        bint optLinearAlgebraInLastBlock
        bint optModifiedLinearAlgebra
        bint optRedTailInLastBlock
        bint optDrawMatrices
        std_string matrixPrefix

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

    int (* pairs_top_sugar)(PBGBStrategy strat)
    PBPolyVector (* someNextDegreeSpolys)(PBGBStrategy strat, int n)
    PBPolyVector (* nextDegreeSpolys)(PBGBStrategy strat)
    PBPolyVector (* small_next_degree_spolys)(PBGBStrategy strat, double f, int n)
    void (* implications)(PBGBStrategy strat, int i)

    PBPoly GB_get_ith_gen "get_ith_gen" (PBGBStrategy strat, int i)

    # non-allocating versions
    PBGBStrategy* PBGBStrategy_construct "Construct<GroebnerStrategy>"(void *mem)
    PBGBStrategy* PBGBStrategy_construct_gbstrategy \
            "Construct_p<GroebnerStrategy, GroebnerStrategy>" \
            (void *mem, PBGBStrategy strat)
    PBPoly* PBPoly_construct_dd \
            "Construct_p<BoolePolynomial, BoolePolyRing::dd_type>" \
            (void *mem, PBDD d)
    void PBGBStrategy_destruct "Destruct<GroebnerStrategy>"(PBGBStrategy *mem)


    PBRedStrategy* PBRedStrategy_construct "Construct<ReductionStrategy>"(void *mem)
    PBRedStrategy* PBRedStrategy_construct_redstrategy \
            "Construct_p<ReductionStrategy, ReductionStrategy>" \
            (void *mem, PBRedStrategy strat)
    void PBRedStrategy_destruct "Destruct<ReductionStrategy>"(PBRedStrategy *mem)

    PBFglmStrategy* PBFglmStrategy_construct "Construct_ppp<FGLMStrategy, BoolePolyRing, BoolePolyRing, PolynomialVector>"(void *mem, PBRing from_ring, PBRing to_ring, PBPolyVector vec)
    void PBFglmStrategy_destruct "Destruct<FGLMStrategy>"(PBFglmStrategy *mem)


   # allocating versions

    PBRedStrategy* PBRedStrategy_new "New<ReductionStrategy>"()
    void PBRedStrategy_delete "Delete<ReductionStrategy>"(PBRedStrategy *mem)

    PBRing get_current_ring "BooleEnv::ring"()

    PBPoly pb_add_up_polynomials "add_up_polynomials" \
            (PBPolyVector v)

    PBPoly pb_nf3 "nf3" (PBRedStrategy strat, PBPoly p, PBMonom m)

    PBPoly pb_red_tail "red_tail" (PBRedStrategy strat, PBPoly p)

    PBPoly pb_interpolate "interpolate" (PBSet z, PBSet o)

    PBPoly pb_interpolate_smallest_lex "interpolate_smallest_lex" \
            (PBSet z, PBSet o)

    PBSet pb_contained_variables_cudd_style "contained_variables_cudd_style" \
            (PBSet m)

    PBSet pb_mod_var_set "mod_var_set" (PBSet a, PBSet v)

    PBPoly pb_mult_fast_sim "mult_fast_sim" (PBPolyVector v)

    #int pb_select1 "select1" (GBStrategy s, PBMonom m)

    PBSet pb_recursively_insert "recursively_insert"(PBNavigator p,
                                                int idx, PBSet m)
    PBPoly pb_ll_red_nf_noredsb "ll_red_nf_noredsb"(PBPoly p,
                                                PBSet reductors)
    PBPoly pb_ll_red_nf_noredsb_single_recursive_call "ll_red_nf_noredsb_single_recursive_call"(PBPoly p,
                                                                                                PBSet reductors)


    PBPoly pb_ll_red_nf "ll_red_nf"(PBPoly p, PBSet reductors)

    PBSet pb_mod_mon_set "mod_mon_set"(PBSet as, PBSet vs)

    PBPolyVector pb_parallel_reduce "parallel_reduce" \
        (PBPolyVector inp, PBGBStrategy strat, int average_steps, double delay_f)

    PBPolyVector pb_gauss_on_polys "gauss_on_polys" (PBPolyVector inp)

    PBPolyVector pb_easy_linear_factors "easy_linear_factors"(PBPoly p)

    PBPoly pb_substitute_variables "substitute_variables<std::vector<BoolePolynomial>, BoolePolynomial>" (PBPolyVector vec, PBPoly poly)

    void pb_append_block "BooleEnv::appendBlock" (int ind)

    void pb_set_variable_name "BooleEnv::setVariableName" \
        (int idx, char *varname)

    char * pb_get_variable_name "BooleEnv::getVariableName" \
        (int idx)

    PBSet pb_random_set "random_set" (PBMonom variables, int length)
    void pb_set_random_seed "set_random_seed" (unsigned int seed)

    #M4RI initialization
    void m4ri_build_all_codes()
    void m4ri_destroy_all_codes()

    void PBPolyVector_set(PBPolyVector v, int i, PBPoly p)
