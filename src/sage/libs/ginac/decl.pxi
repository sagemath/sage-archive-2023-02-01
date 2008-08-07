cdef extern from "ginac_wrap.h":
    ctypedef struct GBasic "basic":
        unsigned int gethash()
        int compare(GBasic other)

    ctypedef struct GSymbol "symbol":
        pass

    GSymbol* GSymbol_construct_str "Construct_p<symbol, char*>" \
            (void *mem, char* m)

    void GSymbol_destruct "Destruct<symbol>"(GSymbol *mem)

    object GSymbol_to_str "_to_PyString<symbol>"(GSymbol *s)

    ctypedef struct GEx "ex":
        unsigned int gethash()
        int compare(GEx other)
        GEx expand(unsigned int opt)
        GEx collect(GEx s, bint dist)

    void GEx_destruct "Destruct<ex>"(GEx *mem)
    GEx* GEx_construct_symbol "Construct_p<ex, symbol>" (void *mem, GSymbol m)
    GEx* GEx_construct_ex "Construct_p<ex, ex>" (void *mem, GEx m)
    GEx* GEx_construct_long "Construct_p<ex, long>" (void *mem, long n)
    GEx* GEx_construct_double "Construct_p<ex, double>" (void *mem, double d)

    object GEx_to_str "_to_PyString<ex>"(GEx *s)

    GEx gadd "ADD_WRAP" (GEx left, GEx right)
    GEx gmul "MUL_WRAP" (GEx left, GEx right)
    GEx gpow "pow" (GEx left, GEx exp)

    GSymbol get_symbol(char* s)
