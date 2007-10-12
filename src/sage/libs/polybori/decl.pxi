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

    ctypedef struct PBDD "struct CDDInterface<CTypes::dd_base>":
        pass

    # non-allocating versions
    void PBDD_destruct "Destruct<CDDInterface<CTypes::dd_base> >"(PBDD *mem)

    ctypedef struct PBNavigator "struct CCuddNavigator":
        PBNavigator (* thenBranch) ()
        PBNavigator (* elseBranch) ()
        int (* value "operator*")()

    # non-allocating versions
    PBNavigator* PBNavigator_construct(void* mem, PBNavigator N)
    void PBNavigator_destruct "Destruct<CCuddNavigator>"(PBNavigator *mem)

    ctypedef struct PBRing "struct BoolePolyRing":
        int (* nVariables)()
        PBDD (* variable)(int n)
        void (*setRingVariableName)(int idx, char *varname)
        void (*activate)()

    ctypedef struct PBSet "struct BooleSet":
        PBNavigator (* navigation)()

    # non-allocating versions
    PBSet* PBSet_construct "Construct<BooleSet>"(void* mem)
    PBSet* PBSet_construct_dd(void* mem, PBDD d)
    void PBSet_destruct "Destruct<BooleSet>"(PBSet *mem)

    ctypedef struct PBMonom "struct BooleMonomial":
        #FIXME need iterator
        int (* deg)()
        int (* hash)()
        comparecodes (* compare)(PBMonom rhs)

    char* PBMonom_to_str "to_str<BooleMonomial>"(PBMonom *p)

    # non-allocating versions
    PBMonom* PBMonom_construct "Construct<BooleMonomial>"(void *mem)
    void PBMonom_destruct "Destruct<BooleMonomial>"(PBMonom *mem)

    ctypedef struct PBPoly "struct BoolePolynomial":
        int (* deg)()
        int (* lmDeg)()
        bint (* isZero)()
        bint (* isOne)()
        bint (* isConstant)()
        PBMonom (* lead)()
        PBDD (* diagram)()
        PBNavigator (* navigation)()

    # allocating versions
    PBRing* PBRing_new "New<BoolePolyRing>"()
    void PBRing_delete "Delete<BoolePolyRing>"(PBRing *mem)

    # non-allocating versions
    PBRing* PBRing_construct (void *mem, int nvars,
                                                        ordercodes order)
    void PBRing_destruct "Destruct<BoolePolyRing>"(PBRing *mem)

    # allocating versions
    PBPoly* PBPoly_new "New<BoolePolynomial>"()
    void PBPoly_delete "Delete<BoolePolynomial>"(PBPoly *mem)

    # non-allocating versions
    PBPoly* PBPoly_construct "Construct<BoolePolynomial>"(void *mem)
    PBPoly PBPoly_construct_dd (void *mem, PBDD d)
    PBPoly PBPoly_construct_pbpoly (void *mem, PBPoly d)
    PBPoly PBPoly_construct_pbmonom (void *mem, PBMonom d)
    PBPoly PBPoly_construct_int (void *mem, int d)
    void PBPoly_destruct "Destruct<BoolePolynomial>"(PBPoly *mem)

    char* PBPoly_to_str "to_str<BoolePolynomial>"(PBPoly *p)

    # PBPoly arithmetic
    PBPoly PBPoly_add (PBPoly left, PBPoly right)
    PBPoly PBPoly_mul (PBPoly left, PBPoly right)

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


    PBSet recursively_insert(PBNavigator p, int idx, PBNavigator m)
    PBPoly ll_red_nf_noredsb(PBPoly p, PBSet reductors)
