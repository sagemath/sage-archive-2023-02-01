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

    # Constants
    GEx g_Pi "Pi"
    GEx g_Catalan "Catalan"
    GEx g_Euler "Euler"

    # Destructor and constructor
    void GEx_destruct "Destruct<ex>"(GEx *mem)
    GEx* GEx_construct_symbol "Construct_p<ex, symbol>" (void *mem, GSymbol m)
    GEx* GEx_construct_ex "Construct_p<ex, ex>" (void *mem, GEx m)
    GEx* GEx_construct_long "Construct_p<ex, long>" (void *mem, long n)
    GEx* GEx_construct_double "Construct_p<ex, double>" (void *mem, double d)

    GEx* GEx_construct_pyobject "ASSIGN_WRAP" (GEx mem, object n)

    object GEx_to_str "_to_PyString<ex>"(GEx *s)

    GEx gadd "ADD_WRAP" (GEx left, GEx right)
    GEx gsub "SUB_WRAP" (GEx left, GEx right)
    GEx gmul "MUL_WRAP" (GEx left, GEx right)
    GEx gdiv "DIV_WRAP" (GEx left, GEx right)
    GEx g_lt "LT_WRAP" (GEx left, GEx right)
    GEx g_eq "EQ_WRAP" (GEx left, GEx right)
    GEx g_gt "GT_WRAP" (GEx left, GEx right)
    GEx g_le "LE_WRAP" (GEx left, GEx right)
    GEx g_ne "NE_WRAP" (GEx left, GEx right)
    GEx g_ge "GE_WRAP" (GEx left, GEx right)
    GEx g_pow "pow" (GEx left, GEx exp)

    GSymbol get_symbol(char* s)

    GEx g_abs "GiNaC::abs" (GEx x)
    GEx g_step "GiNaC::step" (GEx x)	# step function
    GEx g_csgn "GiNaC::csgn" (GEx x)	# complex sign
    GEx g_conjugate "GiNaC::conjugate" (GEx x)	# complex conjugation
    GEx g_real_part "GiNaC::real_part" (GEx x)	# real part
    GEx g_imag_part "GiNaC::imag_part" (GEx x)	# imaginary part
    GEx g_sqrt "GiNaC::sqrt" (GEx x)	# square root (not a GiNaC function, rather an alias for pow(x, numeric(1, 2)))
    GEx g_sin "GiNaC::sin" (GEx x)	# sine
    GEx g_cos "GiNaC::cos" (GEx x)	# cosine
    GEx g_tan "GiNaC::tan" (GEx x)	# tangent
    GEx g_asin "GiNaC::asin" (GEx x)	# inverse sine
    GEx g_acos "GiNaC::acos" (GEx x)	# inverse cosine
    GEx g_atan "GiNaC::atan" (GEx x)	# inverse tangent
    GEx g_atan2 "GiNaC::atan2" (GEx y, GEx x)	# inverse tangent with two arguments
    GEx g_sinh "GiNaC::sinh" (GEx x)	# hyperbolic sine
    GEx g_cosh "GiNaC::cosh" (GEx x)	# hyperbolic cosine
    GEx g_tanh "GiNaC::tanh" (GEx x)	# hyperbolic tangent
    GEx g_asinh "GiNaC::asinh" (GEx x)	# inverse hyperbolic sine
    GEx g_acosh "GiNaC::acosh" (GEx x)	# inverse hyperbolic cosine
    GEx g_atanh "GiNaC::atanh" (GEx x)	# inverse hyperbolic tangent
    GEx g_exp "GiNaC::exp" (GEx x)	# exponential function
    GEx g_log "GiNaC::log" (GEx x)	# natural logarithm
    GEx g_Li2 "GiNaC::Li2" (GEx x)	# dilogarithm
    GEx g_Li "GiNaC::Li" (GEx m, GEx x)	# classical polylogarithm as well as multiple polylogarithm
    GEx g_G "GiNaC::G" (GEx a, GEx y)	# multiple polylogarithm
    GEx g_G2 "GiNaC::G" (GEx a, GEx s, GEx y)	# multiple polylogarithm with explicit signs for the imaginary parts
    GEx g_S "GiNaC::S" (GEx n, GEx p, GEx x)	# Nielsen's generalized polylogarithm
    GEx g_H "GiNaC::H" (GEx m, GEx x)	        # harmonic polylogarithm
    GEx g_zeta "GiNaC::zeta" (GEx m)	        # Riemann's zeta function as well as multiple zeta value
    GEx g_zeta2 "GiNaC::zeta" (GEx m, GEx s)	# alternating Euler sum
    GEx g_zetaderiv "GiNaC::zetaderiv" (GEx n, GEx x)	# derivatives of Riemann's zeta function
    GEx g_tgamma "GiNaC::tgamma" (GEx x)	# gamma function
    GEx g_lgamma "GiNaC::lgamma" (GEx x)	# logarithm of gamma function
    GEx g_beta "GiNaC::beta" (GEx x, GEx y)	# beta function (tgamma*tgamma(y)/tgamma(x+y))
    GEx g_psi "GiNaC::psi" (GEx x)	        # psi (digamma) function
    GEx g_psi2 "GiNaC::psi" (GEx n, GEx x)	# derivatives of psi function (polygamma functions)
    GEx g_factorial "GiNaC::factorial" (GEx n)	# factorial function n!
    GEx g_binomial "GiNaC::binomial" (GEx n, GEx k)	# binomial coefficients
    GEx g_Order "GiNaC::Order" (GEx x)	        # order term function in truncated power series

