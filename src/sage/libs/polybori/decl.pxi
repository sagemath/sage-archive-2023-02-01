cdef extern from "pb_wrap.h":

    #ctypedef size_type "size_type":
    #    pass

    #ctypedef order_type "CTypes::ordercode_type":
    #    pass

    cdef enum ordercodes "COrderEnums::ordercodes":
        lp            "CTypes::lp"
        dlex          "CTypes::dlex"
        dp_asc        "CTypes::dp_asc"
        block_dlex    "CTypes::block_dlex"
        block_dp_asc  "CTypes::block_dp_asc"


    ctypedef struct PBDD "struct CDDInterface<CTypes::dd_base>":
        pass

    # really, this is from NTL/ZZ.h
    ctypedef struct PBRing "struct BoolePolyRing":
        int (* nVariables)()
        PBDD (* variable)(int n)
        void (*setRingVariableName)(int idx, char *varname)

    # Some boiler-plate
    #PBRing* PBRing_new "New<BoolePolyRing>"()
    PBRing* PBRing_construct "BPRing_Construct"(void *mem, int nvars,
                                                        ordercodes order)
    void PBRing_destruct "Destruct<BoolePolyRing>"(PBRing *mem)
    #void PBRing_delete "Delete<BoolePolyRing>"(PBRing *mem)

    ctypedef struct PBPoly "struct BoolePolynomial":
        pass

    #PBPoly* PBPoly_construct "Construct<BoolePolynomial>"(void *mem)
    PBPoly* PBPoly_new "New<BoolePolynomial>"()
    PBPoly* PBPoly_new_dd "BPolyNewDD"(PBDD d)
    PBPoly* PBPoly_new_pbpoly "BPolyNewBPoly"(PBPoly d)
    void PBPoly_delete "Delete<BoolePolynomial>"(PBPoly *mem)

    char* PBPoly_to_str "to_str<BoolePolynomial>"(PBPoly *p)

    # PBPoly arithmetic
    PBPoly* PBPoly_add (PBPoly* left, PBPoly* right)
    PBPoly* PBPoly_mul (PBPoly* left, PBPoly* right)

