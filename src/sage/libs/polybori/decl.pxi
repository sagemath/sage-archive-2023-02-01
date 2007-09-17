cdef extern from "pb_wrap.h":
    cdef enum ordercodes "COrderEnums::ordercodes":
        lp            "CTypes::lp"
        dlex          "CTypes::dlex"
        dp_asc        "CTypes::dp_asc"
        block_dlex    "CTypes::block_dlex"
        block_dp_asc  "CTypes::block_dp_asc"


    ctypedef struct PBDD "struct CDDInterface<CTypes::dd_base>":
        pass

    ctypedef struct PBRing "struct BoolePolyRing":
        int (* nVariables)()
        PBDD (* variable)(int n)
        void (*setRingVariableName)(int idx, char *varname)
        void (*activate)()

    ctypedef struct PBPoly "struct BoolePolynomial":
        pass

    # Some boiler-plate

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

    GBStrategy* GBStrategy_construct "Construct<GroebnerStrategy>"(void *mem)
    void GBStrategy_destruct "Destruct<GroebnerStrategy>"(GBStrategy *mem)
