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

    # non-allocating versions
    void PBDD_destruct "Destruct<CDDInterface<CTypes::dd_base> >"(PBDD *mem)

    ctypedef struct PBRing "BoolePolyRing":
        int (* nVariables)()
        PBDD (* variable)(int n)
        void (*setRingVariableName)(int idx, char *varname)
        void (*activate)()

    ctypedef struct PBSet "BooleSet":
        bint (* emptiness)()
        PBNavigator (* navigation)()

    # non-allocating versions
    PBSet* PBSet_construct "Construct<BooleSet>"(void* mem)
    PBSet* PBSet_construct_dd \
            "Construct_p<BooleSet, BooleSet::dd_type>" (void* mem, PBDD d)
    PBSet* PBSet_construct_pbnav \
            "Construct_p<BooleSet, CCuddNavigator>" (void* mem, PBNavigator d)
    void PBSet_destruct "Destruct<BooleSet>"(PBSet *mem)

    ctypedef struct PBMonomIter "BooleMonomial::const_iterator":
        int (* value "operator*")()
        int (* next "operator++")()
        int (* hash)()

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
    PBMonom* PBMonom_construct_dd \
            "Construct_p<BoolePolynomial, BooleMonomial::dd_type>" (void *mem, PBDD m)
    void PBMonom_destruct "Destruct<BooleMonomial>"(PBMonom *mem)

    ctypedef struct PBPolyIter "BoolePolynomial::ordered_iterator":
        PBMonom (* value "operator*")()
        int (* next "operator++")()
        bint (* equal)(PBPolyIter rhs)

    ctypedef struct PBPoly "BoolePolynomial":
        int (* deg)()
        int (* lmDeg)()
        bint (* isZero)()
        bint (* isOne)()
        bint (* isConstant)()
        bint (* isSingleton)()
        PBMonom (* lead)()
        PBMonom (* usedVariables)()
        PBDD (* diagram)()
        PBNavigator (* navigation)()
        PBPolyIter (* orderedBegin)()
        PBPolyIter (* orderedEnd)()
        void (* iadd "operator+=")(PBPoly right)
        void (* imul "operator*=")(PBPoly right)
        void (* imul_monom "operator*=")(PBMonom right)


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

    ctypedef struct PBPoly_vector "std::vector<BoolePolynomial>":
        int (* size)()
        PBPoly (* get "operator[]")(int)

    ctypedef struct GBStrategy "struct GroebnerStrategy":
        void (* addGeneratorDelayed)(PBPoly)
        void (* symmGB_F2)()
        PBPoly_vector (* minimalize)()

    # non-allocating versions
    GBStrategy* GBStrategy_construct "Construct<GroebnerStrategy>"(void *mem)
    void GBStrategy_destruct "Destruct<GroebnerStrategy>"(GBStrategy *mem)


    PBSet pb_recursively_insert "recursively_insert"(PBNavigator p,
                                                int idx, PBNavigator m)
    PBPoly pb_ll_red_nf_noredsb "ll_red_nf_noredsb"(PBPoly p,
                                                PBSet reductors)
    PBPoly pb_ll_red_nf "ll_red_nf"(PBPoly p, PBSet reductors)
